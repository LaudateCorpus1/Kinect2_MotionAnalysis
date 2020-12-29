# -*- coding: utf-8 -*-
"""
Created on Wed Jul 23 00:49:17 2014

@author: jdl2
"""
import os
import sys
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib import animation

'''LOAD Kinect Tracking Data'''
fileBase    = "Standing_Dual"
file        = "D:\\Kinect2\\Data\\Data10.10.2014\\Victor\\"+fileBase+".bin";
OutFileName = fileBase+".mp4"
Njoints      = 25;
SizeOfSingle = 4;
Nframes      = int(os.path.getsize(file)/(Njoints*3*SizeOfSingle));
Data         = np.fromfile(file,dtype=float32,count=-1)
Data         = Data.reshape((Nframes,3*Njoints))

# First pass through data to get axis bounds
maxVals = np.array([[-Inf, -Inf, -Inf]])
minVals = np.array([[Inf, Inf, Inf]])
for id in range(0,Nframes):
    Frame   = Data[id,:];
    Frame   = Frame.reshape(Njoints,3);
    minVals = np.min(np.concatenate((np.min(Frame,axis=0).reshape(1,3),minVals),axis=0),axis=0).reshape(1,3);
    maxVals = np.max(np.concatenate((np.max(Frame,axis=0).reshape(1,3),maxVals),axis=0),axis=0).reshape(1,3);
    
# Second pass: render animation
fig = plt.figure();
ax1  = fig.add_subplot(1,2,1,projection='3d')
ax1.set_xlim([minVals[0,0],maxVals[0,0]]);
ax1.set_ylim([minVals[0,1],maxVals[0,1]]);
ax1.set_zlim([minVals[0,2],maxVals[0,2]]);

ax2  = fig.add_subplot(1,2,2,projection='3d')
ax2.set_xlim([minVals[0,0],maxVals[0,0]]);
ax2.set_ylim([minVals[0,1],maxVals[0,1]]);
ax2.set_zlim([minVals[0,2],maxVals[0,2]]);

ax1.view_init(-90,90);
ax2.view_init(0,90);
fig.tight_layout();


# Joints as points
Pose1,  = ax1.plot([], [], [], "yo",markersize=5)
# line joining head to tip of left foot
Seg1    = [3,2,20,1,0,12,13,14,15];
Pose1a, = ax1.plot([], [], [], color="green", lw=2)
# line from Spine base to right foot
Seg2    = [0,16,17,18,19];
Pose1b, = ax1.plot([], [], [], color="green", lw=2)
# line from left hand to right hand
Seg3    = [21,7,6,5,4,20,8,9,10,11,23];
Pose1c, = ax1.plot([], [], [], color="green", lw=2)
# left thumb Segment
Seg4    = [22,6];
Pose1d, = ax1.plot([], [], [], color="green", lw=2)
# right thumb Segment
Seg5    = [24,10];
Pose1e, = ax1.plot([], [], [], color="green", lw=2)

Pose2,  = ax2.plot([], [], [], "yo",markersize=5)
Pose2a, = ax2.plot([], [], [], color="green", lw=2)
Pose2b, = ax2.plot([], [], [], color="green", lw=2)
Pose2c, = ax2.plot([], [], [], color="green", lw=2)
Pose2d, = ax2.plot([], [], [], color="green", lw=2)
Pose2e, = ax2.plot([], [], [], color="green", lw=2)

'''
n = 300;
# n = frame counter
Frame = Data[n,:];
Frame = Frame.reshape(Njoints,3);
# update the line data
Pose1.set_data(Frame[:,0], Frame[:,1])
Pose1.set_3d_properties(Frame[:,2])
Pose1a.set_data(Frame[Seg1,0], Frame[Seg1,1])
Pose1a.set_3d_properties(Frame[Seg1,2])
Pose1b.set_data(Frame[Seg2,0], Frame[Seg2,1])
Pose1b.set_3d_properties(Frame[Seg2,2])
Pose1c.set_data(Frame[Seg3,0], Frame[Seg3,1])
Pose1c.set_3d_properties(Frame[Seg3,2])
Pose1d.set_data(Frame[Seg4,0], Frame[Seg4,1])
Pose1d.set_3d_properties(Frame[Seg4,2])
Pose1e.set_data(Frame[Seg5,0], Frame[Seg5,1])
Pose1e.set_3d_properties(Frame[Seg5,2])

Pose2.set_data(Frame[:,0], Frame[:,1])
Pose2.set_3d_properties(Frame[:,2])
Pose2a.set_data(Frame[Seg1,0], Frame[Seg1,1])
Pose2a.set_3d_properties(Frame[Seg1,2])
Pose2b.set_data(Frame[Seg2,0], Frame[Seg2,1])
Pose2b.set_3d_properties(Frame[Seg2,2])
Pose2c.set_data(Frame[Seg3,0], Frame[Seg3,1])
Pose2c.set_3d_properties(Frame[Seg3,2])
Pose2d.set_data(Frame[Seg4,0], Frame[Seg4,1])
Pose2d.set_3d_properties(Frame[Seg4,2])
Pose2e.set_data(Frame[Seg5,0], Frame[Seg5,1])
Pose2e.set_3d_properties(Frame[Seg5,2])

ax1.view_init(-90,90);
ax2.view_init(0,90);
fig.tight_layout();
'''

def init():
    Pose1.set_data([], [])
    Pose1.set_3d_properties([])
    Pose1a.set_data([], [])
    Pose1a.set_3d_properties([])
    Pose1b.set_data([], [])
    Pose1b.set_3d_properties([])
    Pose1c.set_data([], [])
    Pose1c.set_3d_properties([])
    Pose1d.set_data([], [])
    Pose1d.set_3d_properties([])
    Pose1e.set_data([], [])
    Pose1e.set_3d_properties([])
    
    Pose2.set_data([], [])
    Pose2.set_3d_properties([])
    Pose2a.set_data([], [])
    Pose2a.set_3d_properties([])
    Pose2b.set_data([], [])
    Pose2b.set_3d_properties([])
    Pose2c.set_data([], [])
    Pose2c.set_3d_properties([])
    Pose2d.set_data([], [])
    Pose2d.set_3d_properties([])
    Pose2e.set_data([], [])
    Pose2e.set_3d_properties([])
    
def update(n): 
    # n = frame counter
    print "Writing Frame %d out of %d" % (n,Nframes)
    Frame = Data[n,:];
    Frame = Frame.reshape(Njoints,3);
    # update the line data
    Pose1.set_data(Frame[:,0], Frame[:,1])
    Pose1.set_3d_properties(Frame[:,2])
    Pose1a.set_data(Frame[Seg1,0], Frame[Seg1,1])
    Pose1a.set_3d_properties(Frame[Seg1,2])
    Pose1b.set_data(Frame[Seg2,0], Frame[Seg2,1])
    Pose1b.set_3d_properties(Frame[Seg2,2])
    Pose1c.set_data(Frame[Seg3,0], Frame[Seg3,1])
    Pose1c.set_3d_properties(Frame[Seg3,2])
    Pose1d.set_data(Frame[Seg4,0], Frame[Seg4,1])
    Pose1d.set_3d_properties(Frame[Seg4,2])
    Pose1e.set_data(Frame[Seg5,0], Frame[Seg5,1])
    Pose1e.set_3d_properties(Frame[Seg5,2])
    
    Pose2.set_data(Frame[:,0], Frame[:,1])
    Pose2.set_3d_properties(Frame[:,2])
    Pose2a.set_data(Frame[Seg1,0], Frame[Seg1,1])
    Pose2a.set_3d_properties(Frame[Seg1,2])
    Pose2b.set_data(Frame[Seg2,0], Frame[Seg2,1])
    Pose2b.set_3d_properties(Frame[Seg2,2])
    Pose2c.set_data(Frame[Seg3,0], Frame[Seg3,1])
    Pose2c.set_3d_properties(Frame[Seg3,2])
    Pose2d.set_data(Frame[Seg4,0], Frame[Seg4,1])
    Pose2d.set_3d_properties(Frame[Seg4,2])
    Pose2e.set_data(Frame[Seg5,0], Frame[Seg5,1])
    Pose2e.set_3d_properties(Frame[Seg5,2])

print "Writing animaition file..."
sys.stdout.flush()
anim = animation.FuncAnimation(fig, update, init_func=init, frames=Nframes, blit=True)
mywriter     = animation.FFMpegWriter(fps=30,bitrate=2000)


anim.save(OutFileName,writer=mywriter)
print "DONE."
mywriter.finish()

plt.close(fig)