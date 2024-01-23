import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# 重新定义高斯函数
def gaussian(x, y, sigma_x, sigma_y, mu_x, mu_y):
    return np.exp(-((x - mu_x)**2 / (2 * sigma_x**2) + (y - mu_y)**2 / (2 * sigma_y**2)))

# 重新创建网格
x = np.linspace(-5, 5, 100)
y = np.linspace(-5, 5, 100)
X, Y = np.meshgrid(x, y)

# 高斯分布参数
sigma_x, sigma_y = 1.0, 1.0
mu_x, mu_y = 0.0, 0.0

# 计算高斯函数值
Z = gaussian(X, Y, sigma_x, sigma_y, mu_x, mu_y)

# 生成一个LightSource对象
from matplotlib.colors import LightSource
ls = LightSource()

# 使用原有的'hot_r' colormap，但是通过lightsource淡化颜色
new_cmap = ls.blend_hsv('hot_r', fraction=0.3)

# 绘制热力图，使用自定义的colormap
plt.imshow(Z, extent=[-5, 5, -5, 5], origin='lower', norm=LogNorm(), cmap=new_cmap)
plt.colorbar()
plt.title('2D Gaussian Heatmap with Lighter Color')
plt.xlabel('Line rate at t (Mbps)')
plt.ylabel('Line rate at t + d (Mbps)')
plt.savefig('gaussian_heatmap_lighter.png')