import matplotlib.pyplot as plt
import numpy as np

# 示例数据
data = np.random.normal(0, 1, 100)

# 绘制小提琴图
plt.violinplot(data)

# 添加标题和轴标签
plt.title('violin')
plt.xlabel('X')
plt.ylabel('Y')

# 显示图表
plt.show()
plt.savefig('violin.png')  # 保存图表为图片文件
