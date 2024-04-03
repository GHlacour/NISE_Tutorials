def copy_to_multiple_lines(input_file, output_file, num_lines):
    # Open the input file in read mode
    with open(input_file, 'r') as f_input:
        # Read the content of the input file
        content = f_input.readline().strip()

    # Open the output file in write mode
    with open(output_file, 'w') as f_output:
        # Write the content to the output file for the specified number of lines
        for _ in range(num_lines):
            f_output.write(content + '\n')

# Usage
input_file = 'Dipole_1_snapshot.txt'  # Replace with the name of your input file
output_file = 'Dipole.txt'  # Replace with the name of your output file
num_lines = 1000000

copy_to_multiple_lines(input_file, output_file, num_lines)
