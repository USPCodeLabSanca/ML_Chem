import sys
import time

def process_array(array):
    # Simulate some processing
    time.sleep(2)
    return [x * 2 for x in array]

if __name__ == '__main__':
    # Skip the script name (sys.argv[0])
    array = [int(arg) for arg in sys.argv[1:]]
    
    print(f"Process started with array: {array}")
    # Process the array
    result = process_array(array)
    
    print(f"Process completed. Input: {array}, Output: {result}")
