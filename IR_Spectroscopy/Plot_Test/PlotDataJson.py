import json
import matplotlib.pyplot as plt

# Load data from JSON file
with open('test.json', 'r') as file:
    data = json.load(file)

# Extract x and y values from the dataset
wavenumbers = [point['x'] for point in data['datasetColl'][0]['data']]
absorbance = [point['y'] for point in data['datasetColl'][0]['data']]

# Create the plot
plt.figure(figsize=(12, 6))
plt.plot(wavenumbers, absorbance, color='black', linewidth=1.5)
plt.xlabel('Wavenumber (cm⁻¹)')
plt.ylabel('Absorbance')
plt.title('IR Spectrum')

# Remove x-axis inversion
# plt.gca().invert_xaxis()

# Invert y-axis for conventional IR representation
plt.gca().invert_yaxis()

# Set y-axis limits based on calibration points and actual data
calibration_points = data['axesColl'][0]['calibrationPoints']
y_min_cal = min(float(point['dy']) for point in calibration_points)
y_max_cal = max(float(point['dy']) for point in calibration_points)
y_min_data = min(absorbance)
y_max_data = max(absorbance)

y_min = min(y_min_cal, y_min_data)
y_max = max(y_max_cal, y_max_data)

# Add some padding to y-axis
y_range = y_max - y_min
y_padding = 0.05 * y_range
y_min -= y_padding
y_max += y_padding

# Set x-axis limits based on actual data points
x_min = min(wavenumbers)
x_max = max(wavenumbers)

plt.xlim(x_min, x_max)  # Ascending x-axis
plt.ylim(y_max, y_min)  # Inverted y-axis

# Display the plot
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()

# Adjust the plot to be centered
plt.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.1)

plt.show()
