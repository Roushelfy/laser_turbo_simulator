import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# 定义高斯函数
def gaussian(x, y, sigma_x, sigma_y, mu_x, mu_y):
    return np.exp(-((x+ - mu_x)**2 / (2 * sigma_x**2) + (y - mu_y)**2 / (2 * sigma_y**2)))

# 创建网格
x = np.linspace(-5, 5, 100)
y = np.linspace(-5, 5, 100)
X, Y = np.meshgrid(x, y)

# 高斯分布参数
sigma_x, sigma_y = 1.5, 1.5
mu_x, mu_y = 1.5, 1.5

# 计算具有噪声的高斯函数值
Z = gaussian(X, Y, sigma_x, sigma_y, mu_x, mu_y)*(1+0.1*np.random.randn(100,100))

# 绘制热力图，使用LogNorm来得到对数刻度的色彩映射
plt.imshow(Z, extent=[-5, 5, -5, 5], origin='lower', norm=LogNorm(), cmap='Reds')
plt.colorbar()
plt.title('2D Gaussian Heatmap with Logarithmic Color Scale')

# 画虚线格子
plt.grid(linestyle='--')
plt.show()
plt.savefig('gaussian_heatmap.pdf')
