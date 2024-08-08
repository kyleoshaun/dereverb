function plot_LIVELab(idx)

load('spklocs.mat') % Load locations of LIVELab speakers; given in units of metre, relative to the bottom back-left corner of the LIVELab

figure

% Plot absolute speaker locations
plot3(spklocs(:,1),spklocs(:,2),spklocs(:,3),'o','markersize',14)
hold on
for lp=1:length(spknmbrs)
    if isempty(find(idx==lp)) == false
        text(spklocs(lp,1),spklocs(lp,2),spklocs(lp,3),spknmbrs{lp},'fontsize',7,'horizontalalignment','center',  'Color','red')
    else
        text(spklocs(lp,1),spklocs(lp,2),spklocs(lp,3),spknmbrs{lp},'fontsize',7,'horizontalalignment','center')
    end
end

ylim([0 15.8])
ylabel('y (m) - back wall to front wall')
xlim([0 11.35])
xlabel('x (m) - left wall to right wall')
zlim([0 6])
zlabel('z (m) - floor to ceiling')
title('Absolute loudspeaker and microphone locations')
grid on

load('micloc_example.mat') % Load location of microphone or KEMAR (vertex) for specific recording

plot3(micloc(1),micloc(2),micloc(3),'r.','markersize',6) % Plot microphone or KEMAR vertex location


spklocs_rel = spklocs - micloc; % calculate loudspeaker positions relative to the microphone

[az,el,r] = cart2sph(spklocs_rel(:,2),spklocs_rel(:,1),spklocs_rel(:,3)); % convert relative positions to spherical coordinates,
                                                                          % switching the x and y axes so that 0 degrees is the
                                                                          % forward direction, with positive angles to the right

azd = az/pi*180; % convert azimuth values to degrees
eld = el/pi*180; % convert elevation values to degrees

figure

% Plot speaker locations relative to microphone
plot3(spklocs_rel(:,1),spklocs_rel(:,2),spklocs_rel(:,3),'o','markersize',14)
hold on
for lp=1:length(spknmbrs)
    if isempty(find(idx==lp)) == false
        text(spklocs_rel(lp,1),spklocs_rel(lp,2),spklocs_rel(lp,3),spknmbrs{lp},'fontsize',7,'horizontalalignment','center', 'Color','red')
    else
        text(spklocs_rel(lp,1),spklocs_rel(lp,2),spklocs_rel(lp,3),spknmbrs{lp},'fontsize',7,'horizontalalignment','center')
    end
end

plot3(0,0,0,'r.','markersize',6) % Plot microphone or KEMAR vertex location

ylabel('y\_rel (m) - back wall to front wall')
xlabel('x\_rel (m) - left wall to right wall')
zlabel('z\_rel (m) - floor to ceiling')
title('Loudspeaker locations relative to microphone')
grid on

spknmbrs_mat = cell2mat(spknmbrs'); % convert speaker number cell array to a matrix 

sphcordmat = [str2num(spknmbrs_mat(:,2:end)) azd eld r];

save sphcord_mat sphcordmat

end