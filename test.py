from __init__ import *

arry = empty((3, 10, 10))
for k in range(3):
    
f, ax = plt.subplots(figsize=(6, 6))

x = rand(100)
y = rand(100)
sns.kdeplot(x, y, ax=ax, cmap="Blues")
x = rand(100)
y = rand(100)
sns.kdeplot(x, y, ax=ax, cmap="Greens")
x = rand(100)
y = rand(100)
sns.kdeplot(x, y, ax=ax, cmap="Reds")
#sns.rugplot(x, color="g", ax=ax)
#sns.rugplot(y, vertical=True, ax=ax);

plt.savefig('test.png')
plt.show()

