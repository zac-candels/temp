import os
import re
import subprocess

home_dir = os.getcwd()  # assuming you're running this from your home directory
output_file = os.path.join(home_dir, "velLandscape.txt")

with open(output_file, "w") as out:
    for dirname in sorted(os.listdir(home_dir)):
        if not dirname.startswith("test_CA"):
            continue
        
        # Extract x and y from directory name
        match = re.match(r"test_CA(\d+)_postFrac([\d.]+)", dirname)
        if not match:
            continue
        x, y = match.groups()
        
        script_path = os.path.join(home_dir, dirname, "centroid_velocity.py")
        if not os.path.isfile(script_path):
            print(f"⚠️ Skipping {dirname}, no centroid_velocity.py found")
            continue
        
        # Run the python script and capture output
        try:
            result = subprocess.run(
                ["python3", script_path],
                cwd=os.path.join(home_dir, dirname),
                capture_output=True,
                text=True,
                check=True
            )
            # Find the line containing vel

            #print(f"--- Output from {dirname} ---\n{result.stdout}\n--- STDERR ---\n{result.stderr}")

            vel_match = re.search(r"vel\s*=\s*([-\d.eE]+)", result.stdout)
            if vel_match:
                vel = vel_match.group(1)
                out.write(f"{x},{y},{float(vel):.16e}\n")
                print(f"✅ {dirname}: vel={vel}")
            else:
                print(f"⚠️ No velocity found in {dirname}")
        except subprocess.CalledProcessError as e:
            print(f"❌ Error running {dirname}: {e}")
