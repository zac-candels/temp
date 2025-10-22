import os
import shutil
import sys

# -----------------------------
# Configure your paths here
# -----------------------------
src_dir = "/work/sc139/sc139/s2767572/centroidVelLandscape_lx128_pts90"
dst_dir = "/work/sc139/sc139/s2767572/centroidVelLandscape_lx256_pts90"  # somewhere you can write

# File extensions to exclude
exclude_ext = {".vtk", ".dat", ".log", ".mat", ".mp4", ".out"}

# -----------------------------
# Copy function
# -----------------------------
def copy_directory(src, dst, exclude_ext=None):
    exclude_ext = exclude_ext or set()
    
    if not os.path.exists(src):
        raise FileNotFoundError(f"Source directory '{src}' does not exist")
    
    os.makedirs(dst, exist_ok=True)

    for root, dirs, files in os.walk(src):
        rel_path = os.path.relpath(root, src)
        dest_root = os.path.join(dst, rel_path)
        os.makedirs(dest_root, exist_ok=True)

        for file in files:
            if not any(file.endswith(ext) for ext in exclude_ext):
                src_file = os.path.join(root, file)
                dst_file = os.path.join(dest_root, file)
                try:
                    shutil.copy2(src_file, dst_file)
                except PermissionError:
                    print(f"Skipping '{dst_file}': Permission denied", file=sys.stderr)
                except Exception as e:
                    print(f"Error copying '{dst_file}': {e}", file=sys.stderr)


# -----------------------------
# Run the copy
# -----------------------------
copy_directory(src_dir, dst_dir, exclude_ext)
print(f"Finished copying from '{src_dir}' to '{dst_dir}' (excluded: {exclude_ext})")
