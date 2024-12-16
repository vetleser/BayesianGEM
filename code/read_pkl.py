import pickle
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')


# Specify the path to your .pkl file
file_path = "../results/crowdingDE/smcevo_gem_1.0_0.999_0.pkl"

# Load the file
with open(file_path, "rb") as file:
    data = pickle.load(file)

# Inspect the data
print(type(data))  # Check the type of the data
print(data)        # Print the content (if it's not too large)
