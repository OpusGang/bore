#!/usr/bin/env python3

import sys
import shutil
import os

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {os.path.basename(sys.argv[0])} <source_file_to_install> <destination_directory>", file=sys.stderr)
        sys.exit(1)

    source_file_path = sys.argv[1]  # Meson will substitute the path to the built libdescale.so/dll here
    destination_dir = sys.argv[2]   # The target installation directory (e.g., /usr/lib/vapoursynth)

    if not os.path.exists(source_file_path):
        print(f"Error: Source file '{source_file_path}' does not exist at install time. This should not happen if the build succeeded.", file=sys.stderr)
        sys.exit(1)

    # Ensure the destination directory exists
    if not os.path.isdir(destination_dir):
        try:
            # exist_ok=True: don't raise an error if the directory already exists
            os.makedirs(destination_dir, exist_ok=True)
            print(f"Created directory: {destination_dir}")
        except OSError as e:
            print(f"Error: Could not create destination directory '{destination_dir}': {e}", file=sys.stderr)
            sys.exit(1)

    destination_file_name = os.path.basename(source_file_path)
    destination_file_path = os.path.join(destination_dir, destination_file_name)

    try:
        shutil.copy2(source_file_path, destination_file_path) # copy2 preserves metadata like permissions
        print(f"Installed '{source_file_path}' to '{destination_file_path}'")
    except Exception as e:
        print(f"Error: Could not copy '{source_file_path}' to '{destination_file_path}': {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
