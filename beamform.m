% Perform a basic diverging wave beamforming on the complete data set
function [rf_focused] = beamform(rf_fsa,r,rx_pos,x,z)

% Pre-calculate distances
[xg,zg]=meshgrid(x,z);
xg=xg(:);
zg=zg(:);
delays=zeros(length(z)*length(x),size(rx_pos,1),'single');
for i=1:size(rx_pos,1)
    delays(:,i)=sqrt((xg-rx_pos(i,1)).^2+(zg-rx_pos(i,3)).^2);
end

%Loop over element pairs and focus, sum
rf_focused=zeros(length(z)*length(x),1,'single');
for i=1:size(delays,2)
    for j=1:size(delays,2)
        rf_focused=rf_focused+interp1(r,rf_fsa(:,i,j),delays(:,i)+delays(:,j));
    end
end
rf_focused=reshape(rf_focused,length(z),length(x));