#!/usr/bin/env python3

import os
import shutil
import sys
from pathlib import Path


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {Path(sys.argv[0]).name} <source_file_to_install> <destination_directory>", file=sys.stderr)
        sys.exit(1)

    specified_destdir = os.environ.get("DESTDIR", None)
    source_file_path = Path(sys.argv[1])  # Meson will substitute the path to the built libdescale.so/dll here
    destination_dir = Path(sys.argv[2])  # The target installation directory (e.g., /usr/lib/vapoursynth)

    if specified_destdir:
        specified_destdir = Path(specified_destdir)
        destination_dir = specified_destdir / destination_dir.relative_to(destination_dir.anchor)

    if not source_file_path.exists():
        print(
            f"Error: Source file '{source_file_path}' does not exist at install time. This should not happen if the build succeeded.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Ensure the destination directory exists
    if not destination_dir.is_dir():
        try:
            # exist_ok=True: don't raise an error if the directory already exists
            destination_dir.mkdir(parents=True, exist_ok=True)
            print(f"Created directory: {destination_dir}")
        except OSError as e:
            print(f"Error: Could not create destination directory '{destination_dir}': {e}", file=sys.stderr)
            sys.exit(1)

    destination_file_name = source_file_path.name
    destination_file_path = destination_dir / destination_file_name

    try:
        shutil.copy2(source_file_path, destination_file_path)  # copy2 preserves metadata like permissions
        print(f"Installed '{source_file_path}' to '{destination_file_path}'")
    except Exception as e:
        print(f"Error: Could not copy '{source_file_path}' to '{destination_file_path}': {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
