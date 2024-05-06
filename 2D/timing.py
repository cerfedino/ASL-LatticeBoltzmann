#This file runs the main.cpp file multiple times and records 
#When testing, don't forget to update the fields marked with "CHANGE"
#the result will be saved in the results folder
import subprocess

# Number of times to run the C++ code
num_runs = 1

#different input sets to be explored for x, y resolution
inputs = ["8 2", "16 4", "32 8", "64 16", "128 32", "256 64", "400 100"]

#CHANGE FROM HERE
name = "eb4e798a19b616a8f674b32f849117cc935f9365"
description = "Misc optmization"
numb_flops = 10000
#CHANGE TO HERE

def run_main_with_input(cpp_input):
    

    # List to store the outputs
    outputs = []

    # Run the C++ code num_runs times
    for run_numb in range(num_runs):

        # Run the compiled C++ program and capture its output
        result = subprocess.run(["./main.exe"],input = cpp_input.encode(), stdout=subprocess.PIPE)

        #extract the number of cycles
        print(result.stdout.decode().strip())
        if(len(result.stdout.decode().strip().split()) < 2 ): return -1
        output = result.stdout.decode().strip().split()[2]
        print("Curent run:", run_numb, "out of", num_runs, "with that many cycles:", output)
        outputs.append(int(output))

    # Calculate the average
    average_output = sum(outputs) / num_runs

    print("Average output:", average_output)
    return average_output

def write_data_to_file(result, filename):
    try:
        # Open the file in write mode (creates the file if it doesn't exist)
        with open(filename, 'a') as file:
            # Write the result to the file
            file.write(str(result))
            file.write('\n')
        print("Result successfully saved to", filename)
    except IOError:
        # Handle file IO errors
        print("Error: Unable to write to file", filename)


def clear_file(filename):
    with open(filename, 'w') as file:
            file.write(" ")
    
print(description)

results = []
pos = 0
for i in inputs:
    print("Current input: ", i)
    results.append(run_main_with_input(i))

result_filename = "results/" + name + "_result.txt"

clear_file(result_filename)
write_data_to_file("Name: "+name, result_filename)
write_data_to_file("Description: "+description, result_filename)
write_data_to_file("Inputs: ", result_filename)
write_data_to_file(inputs, result_filename)
write_data_to_file("Results: ", result_filename)
write_data_to_file(results, result_filename)
write_data_to_file("Average performance with flops as "+str(numb_flops), result_filename)
performance = [element / numb_flops for element in results]

write_data_to_file(performance, result_filename)

print("Done, result is saved in", result_filename)

