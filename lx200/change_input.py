import os
import re

# Loop over all directories in the current working directory
for dirname in os.listdir("."):
    if not dirname.startswith("test_CA"):
        continue
    
    # Extract x and y from the directory name
    # Format: test_CAx_postFracy
    match = re.match(r"test_CA(\d+)_postFrac([\d.]+)", dirname)
    if not match:
        continue
    
    x, y = match.groups()
    
    input_file = os.path.join(dirname, "input.txt")
    if not os.path.isfile(input_file):
        print(f"⚠️ Skipping {dirname}, no input.txt found")
        continue
    
    # Read file
    with open(input_file, "r") as f:
        lines = f.readlines()
    
    # Modify lines
    new_lines = []
    for line in lines:
        if line.startswith("theta="):
            new_lines.append(f"theta={x} #contact angle\n")
        elif line.startswith("postfraction="):
            new_lines.append(f"postfraction={y} #number of posts in the x direction\n")
        elif line.startswith('datadir="./"'):
            new_lines.append('datadir="./data" #data directory\n')
        elif line.startswith('timesteps='):
            new_lines.append('timesteps=10000000 \n')
        elif line.startswith('saveInterval='):
            new_lines.append('saveInterval=10000 \n')
        elif line.startswith('bodyforcex='):
            new_lines.append('bodyforcex=0.0000006 \n')
        elif line.startswith('lx='):
            new_lines.append('lx=250 \n')
        elif line.startswith('ly='):
            new_lines.append('ly=150 \n')            
        elif line.startswith('posx='):
            new_lines.append('posx=125\n')
        elif line.startswith('radius='):
            new_lines.append('radius=55\n')
        else:
            new_lines.append(line)
    
    # Write back
    with open(input_file, "w") as f:
        f.writelines(new_lines)
    
    print(f"✅ Updated {input_file} with theta={x}, postfraction={y}, datadir=./data")
