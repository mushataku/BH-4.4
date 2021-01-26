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
beta_deg = 15.0
beta = beta_deg/180.0*pi

Rmin,Rmax,dr = 12.0, 15.0, 0.2
r_list = np.arange(Rmin,Rmax+dr,dr)
N_ring = len(r_list)
# 各 ring に含まれる粒子数のリスト
Num_list = [int(50*(r_list[i]/10)) for i in range(N_ring)]

skip = N_ring//5
R_Tg = 10.0
Tg = 2*pi*R_Tg**1.5
###########################

################### variable ##################
SAVE_PATH = './figs/ani_a%.1fbeta%.1f_r%.1f-r%.1f'%(a_spin, beta_deg, Rmin, Rmax)
DIR_PATH = "./output/a%.1fbeta%.1f/a%.1fbeta%.1f"%(a_spin, beta_deg, a_spin, beta_deg)
T = pd.read_csv("%sr%.1f.csv"%(DIR_PATH, Rmin))["time"]
################### variable ##################


#########################描画のための関数#########################

# frame 枚目の写真を作るデータを取得
def get_data():
  ring_list = [] #ring_list[ring][particle][xyz][frame]
  for i_ring in range(N_ring):
    df = pd.read_csv("%sr%.1f.csv"%(DIR_PATH, r_list[i_ring]))
    xyz_list = []
    for i_particle in range(Num_list[i_ring]):
      x = df["x%d"%i_particle]
      y = df["y%d"%i_particle]
      z = df["z%d"%i_particle]
      xyz_list.append([x,y,z])
    ring_list.append(xyz_list)
  return ring_list

# リング i_ring の各運動量ベクトルを返す
def get_L():
  L_list = [] # L_list[ring][xyz][frame]
  for i_ring in range(N_ring):
    if(i_ring%skip != skip-1): continue
    # その半径の ring の 0 番目の粒子の角運動量ベクトルを求める
    df = pd.read_csv("%sr%.1f.csv"%(DIR_PATH, r_list[i_ring]))
    x = 3*np.array(df["x%d"%0])
    y = 3*np.array(df["y%d"%0])
    z = 3*np.array(df["z%d"%0])
    vx = 3*np.array(df["vx%d"%0])
    vy = 3*np.array(df["vy%d"%0])
    vz = 3*np.array(df["vz%d"%0])
    L_list.append([y*vz-z*vy, z*vx-x*vz, x*vy-y*vx])
  return L_list

# 初期画像を設定
def init_ax(ax,ring_list):

  ax.set(xlim=(-Rmax,Rmax),ylim=(-Rmax,Rmax),zlim=(-Rmax/3,Rmax/3))
  # ax.set(xlim=(-Rmax,Rmax),ylim=(-Rmax,Rmax),zlim=(-Rmax,Rmax))

  # ring_list[ring][particle][xyz][frame]
  line_list = []
  for i_ring, R in enumerate(r_list):
    R_norm = (R - Rmin)/(Rmax-Rmin)
    for i_particle in range(Num_list[i_ring]):
      x = ring_list[i_ring][i_particle][0][0]
      y = ring_list[i_ring][i_particle][1][0]
      z = ring_list[i_ring][i_particle][2][0]
      if(i_particle == 0):
        line, = ax.plot([x],[y],[z],'o',color=cm.jet(R_norm), alpha=1.0 ,label="r = %.1f"%R)
      else:
        line, = ax.plot([x],[y],[z],'o', color=cm.jet(R_norm), alpha=0.5)
      line_list.append(line)
  
  ax.legend()
  
  # xy 平面と z 軸のプロット
  x = np.linspace(-Rmax, Rmax, 5)
  y = np.linspace(-Rmax, Rmax, 5)
  X, Y = np.meshgrid(x, y)
  Z = 0*X
  ax.plot_surface(X, Y, Z, alpha=0.3)
  ax.plot([0,0,0],[0,0,0],[-Rmax,0,Rmax],"-", color="black")
  return line_list

# 法線ベクトル
def L_initial(ax, L_list):
  line_L_list = []
  for i,L in enumerate(L_list):
    # L_list[ring][xyz][frame]
    R_norm = (r_list[i] - Rmin)/(Rmax-Rmin)
    line, = ax.plot([-L[0][0],0,L[0][0]],[-L[1][0],0,L[1][0]],[-L[2][0],0,L[2][0]],'-', color=cm.jet(R_norm))
    line_L_list.append(line)

  return line_L_list

# データ更新
def reset_data(line_list,line_L_list,frame,ring_list):
  cnt = 0
  # 位置の更新
  for i_ring in range(N_ring):
    for i_particle in range(Num_list[i_ring]):
      line_list[cnt].set_data(ring_list[i_ring][i_particle][0][frame],ring_list[i_ring][i_particle][1][frame])
      line_list[cnt].set_3d_properties(ring_list[i_ring][i_particle][2][frame])
      cnt += 1
  
  # 角運動量ベクトルの更新
  for i,L in enumerate(L_list):
    line_L_list[i].set_data([-L[0][frame],0,L[0][frame]],[-L[1][frame],0,L[1][frame]])
    line_L_list[i].set_3d_properties([-L[2][frame],0,L[2][frame]])

  

#########################描画のための関数#########################



################################################################
########################### main ###############################
################################################################

# Attaching 3D axis to the figure
fig = plt.figure(figsize=(8,6))
ax = fig.gca(projection='3d')

time_text = fig.text(0.01, 0.99, '', size=20, color="white", horizontalalignment='left',
            verticalalignment='top', backgroundcolor='black')
fig.text(0, 0.01, r"a = %.2f"%a_spin,
          backgroundcolor="black",color="white", size=20)

#### アニメの初期画像生成

ring_list = get_data() #ring_list[ring][particle][xyz][frame]
line_list = init_ax(ax,ring_list)
L_list = get_L()
line_L_list = L_initial(ax,L_list)
time_text.set_text("time = 0.000")

#### 画像更新用関数
def animate(frame):
  if(frame == 0):
    return
  print("frame:",frame)
  reset_data(line_list,  line_L_list, frame, ring_list)
  time_text.set_text(r"time = %.3e Tg for r=%.1f"%(T[frame]/Tg, R_Tg))

ani = FuncAnimation(fig, animate, frames=96
              , interval=200, repeat=True, blit=False)


if(SAVE == 0):
    plt.show()
if(SAVE == 1):
    ani.save(SAVE_PATH+".mp4", writer="ffmpeg", fps=5)
if(SAVE == 2):
    ani.save(SAVE_PATH+".gif", writer='pillow')


