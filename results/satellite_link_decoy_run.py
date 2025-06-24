# -*- coding: utf-8 -*-


import numpy as np
import satellite_link_decoy
import satellite


sat = satellite.Satellite(ra=200e6, time=1)
sat.setup_channels(0)

channel_downlink = sat.downlink_channel

loss_avg = np.load('pass_loss_down.npy')
loss_std = np.load('pass_loss_down_std.npy')


sat_link = satellite_link_decoy.SatelliteLinkDecoy(channel_downlink)

print(sat_link.run_pnt(10, 1))
