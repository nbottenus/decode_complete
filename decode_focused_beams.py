import numpy as np


def decode(rf, delays):
    """
    Decode RF channel data from a set of focused transmit beams into the
    complete data set using the applied focal timings
    rf_fsa = decode(rf,delays);

    Parameters:
    rf - RF data (time sample x receive channel x transmit event)
    transmit_delays - Transmit focal delays in samples (transmit event x
    transmit element)

    Returns:
    rf_fsa - Decoded complete data set (time sample x receive channel x
    transmit element)

    Author: Nick Bottenus
    Contact: nick.bottenus@duke.edu
    """

    RF = np.fft.rfft(rf, axis=0)
    n_samples = RF.shape[0]
    FSA = np.empty([n_samples, RF.shape[1], RF.shape[1]], dtype='complex64')
    Hinv = np.empty(delays.shape, dtype='complex64')
    for frequency in range(0, n_samples):
        phase = 2*np.pi*frequency/n_samples/2*delays
        # Hinv = np.exp(1j*phase)
        Hinv = np.cos(phase)+1j*np.sin(phase)
        FSA[frequency, :, :] = np.dot(RF[frequency, :, :], Hinv)
    fsa = np.fft.irfft(FSA, axis=0)
    return fsa


def beamform(fsa, r, rx_pos, x, z):
    xg, zg = np.meshgrid(x, z)
    xg = xg.flatten()
    zg = zg.flatten()
    delays = np.empty([len(xg), rx_pos.shape[0]], dtype='float32')
    for i, pos in enumerate(rx_pos):
        delays[:, i] = np.sqrt(np.square(xg - pos[0]) +
                               np.square(pos[1]) +
                               np.square(zg - pos[2]))
    r = r.squeeze()
    rf_focused = np.zeros([len(xg)], dtype='float32')
    for i_rx, rx in enumerate(delays.T):
        for i_tx, tx in enumerate(delays.T):
            rf_focused = rf_focused + np.interp(tx + rx, r, fsa[:, i_rx, i_tx])
    rf_focused = rf_focused.reshape([len(z), len(x)])
    return rf_focused


def main():
    import scipy.io as sio
    import scipy.signal as ssi
    import matplotlib.pyplot as plt

    # Load data
    data = sio.loadmat('sample_data.mat')

    # Optionally reduce the number of transmit events in the saved data
    # ds = 6
    # data['rf']=data['rf'][:,:,0::ds]
    # data['transmit_delays'] = data['transmit_delays'][0::ds,:]

    # Perform decoding
    fsa = decode(data['rf'], data['transmit_delays'])

    # Unpack data parameters
    rx_pos = data['params'][0, 0]['rx_pos']
    fs = data['params'][0, 0]['fs']
    t0 = data['params'][0, 0]['t0']
    c = data['params'][0, 0]['c']
    r = (t0+np.array(range(0, fsa.shape[0])))/fs*c

    # Beamforming parameters
    x = np.linspace(-15, 15, num=200, endpoint=True)/1000
    z = np.linspace(15, 60, num=500, endpoint=True)/1000

    # Perform beamforming
    rf_focused = beamform(fsa, r, rx_pos, x, z)
    b, a = ssi.butter(2, 500e3/(fs/2), btype='high')
    rf_focused = ssi.lfilter(b, a, rf_focused, axis=0)

    # Show the output
    env = np.abs(ssi.hilbert(rf_focused, axis=0))
    env = env/env.max()
    plt.imshow(20*np.log10(env), cmap='gray',
               extent=[x[0]*1e3, x[-1]*1e3, z[-1]*1e3, z[0]*1e3],
               vmin=-50, vmax=0)
    plt.xlabel('Lateral (mm)')
    plt.ylabel('Axial (mm)')
    plt.title('Recovered complete data set')
    plt.show()

if __name__ == '__main__':
    main()

