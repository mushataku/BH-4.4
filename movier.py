# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from math import sqrt, pi
from matplotlib.animation import FuncAnimation
from numpy import exp, pi, sin, cos, cosh, abs
from matplotlib.ticker import ScalarFormatter
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D 
import os

##########CONFIG###########
# 0:plt 1:mp4 2:gif
SAVE = 1

a_spin = 1.0
N_particle = 50
N_color = N_particle
r0 = 10.0
Tg = 2*pi*r0**1.5
###########################

################### variable ##################
FILE_PATH = './figs/animation'
df = pd.read_csv("./output/position.csv")
T = df["time"]
################### variable ##################


#########################描画のための関数#########################

# frame 枚目の写真を作るデータを取得
def get_data(df):
  xyz_list = []
  for i in range(N_particle):
    x = df["x%d"%i]
    y = df["y%d"%i]
    z = df["z%d"%i]
    xyz_list.append([x,y,z])
  return xyz_list

# 初期画像を設定
def init_ax(ax,xyz_list):

  ax.set(xlim=(-r0,r0),ylim=(-r0,r0),zlim=(-2,2))
  ax.set(xlim=(-r0,r0),ylim=(-r0,r0),zlim=(-r0,r0))

  line_list = []
  for i in range(N_particle):
    x = xyz_list[i][0][0]
    y = xyz_list[i][1][0]
    z = xyz_list[i][2][0]
    if(i%N_color == 0):
      line, = ax.plot([x],[y],[z],'o',color='r')
    else:
      line, = ax.plot([x],[y],[z],'o', color='black')
    line_list.append(line)
  # ax.legend()

  x = np.linspace(-r0, r0, 5)
  y = np.linspace(-r0, r0, 5)
  X, Y = np.meshgrid(x, y)
  Z = 0*X
  ax.plot_surface(X, Y, Z, alpha=0.3)
  ax.plot([0,0,0],[0,0,0],[-r0,0,r0],"-", color="black")
  return line_list

def normal_ax(ax,xyz_list):
  x1 = xyz_list[0][0][0]
  y1 = xyz_list[0][1][0]
  z1 = xyz_list[0][2][0]
  x2 = xyz_list[1][0][0]
  y2 = xyz_list[1][1][0]
  z2 = xyz_list[1][2][0]

  x = y1*z2-y2*z1
  y = z1*x2-z2*x1
  z = x1*y2-x2*y1
  if(z < 0):
    x = -x
    y = -y
    z = -z
  line, = ax.plot([-x,0,x],[-y,0,y],[-z,0,z],'-', color='black')
  return line

# データ更新
def reset_data(line_list,line_n,frame,xyz_list):
  for i in range(len(line_list)):
    line_list[i].set_data(xyz_list[i][0][frame],xyz_list[i][1][frame])
    line_list[i].set_3d_properties(xyz_list[i][2][frame])
  
  x1 = xyz_list[0][0][frame]
  y1 = xyz_list[0][1][frame]
  z1 = xyz_list[0][2][frame]
  x2 = xyz_list[1][0][frame]
  y2 = xyz_list[1][1][frame]
  z2 = xyz_list[1][2][frame]

  x = (y1*z2-y2*z1)
  y = (z1*x2-z2*x1)
  z = (x1*y2-x2*y1)
  if(z < 0):
    x = -x
    y = -y
    z = -z
  
  

  line_n.set_data([-x,0,x],[-y,0,y])
  line_n.set_3d_properties([-z,0,z])  

#########################描画のための関数#########################



################################################################
########################### main ###############################
################################################################

# Attaching 3D axis to the figure
fig = plt.figure(figsize=(8,6))
# fig.subplots_adjust(left=0.2)
ax = fig.gca(projection='3d')

time_text = fig.text(0.01, 0.99, '', size=20, color="white", horizontalalignment='left',
            verticalalignment='top', backgroundcolor='black')
fig.text(0, 0.01, r"a = %.2f"%a_spin,
          backgroundcolor="black",color="white", size=20)

#### アニメの初期画像生成
xyz_list = get_data(df)
line_list = init_ax(ax,xyz_list)
line_n = normal_ax(ax,xyz_list)
time_text.set_text("time = 0.000")

#### 画像更新用関数
def animate(frame):
  if(frame == 0):
    return
  print("frame:",frame)
  reset_data(line_list,  line_n, frame, xyz_list)
  time_text.set_text(r"time = %.3e Tg"%(T[frame]/Tg))

ani = FuncAnimation(fig, animate, frames=96
              , interval=200, repeat=True, blit=False)


if(SAVE == 0):
    plt.show()
if(SAVE == 1):
    ani.save(FILE_PATH+".mp4", writer="ffmpeg", fps=5)
if(SAVE == 2):
    ani.save(FILE_PATH+".gif", writer='pillow')


