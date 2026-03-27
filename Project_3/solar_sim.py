import pandas as pd
import matplotlib.pyplot as plt


data = pd.read_csv('orbit_data.csv')
sun = data[data['object'] == 'Sun']
jupiter = data[data['object'] == 'Jupiter']
earth = data[data['object'] == 'Earth']


plt.scatter(sun['x'], sun['y'], color = "red", s= 50, label = "Sun" )

plt.plot(earth['x'], earth['y'], color = "green" , label = "Earth")
plt.plot(jupiter['x'], jupiter['y'], color = "brown", label = "Jupiter")


plt.xlabel('x (AU)')
plt.ylabel('y (AU)')
plt.title('Three Body Solar System Simulation')
plt.axis('equal')
plt.grid(True)
plt.legend()
plt.show()