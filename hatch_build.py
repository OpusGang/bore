import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from hatchling.builders.hooks.plugin.interface import BuildHookInterface
from packaging import tags


class CustomHook(BuildHookInterface[Any]):
    source_dir = Path("build")
    target_dir = Path("vapoursynth/plugins")

    def initialize(self, version: str, build_data: dict[str, Any]) -> None:
        build_data["pure_python"] = False

        # the plugin doesn't depend on the Python version so tag py3-none-<platform>
        # instead of the default cpXY-cpXY-<platform> to ship one wheel per platform
        build_data["tag"] = f"py3-none-{next(tags.platform_tags())}"

        # find VapourSynth headers via the pkg-config file shipped in the vapoursynth wheel
        env = os.environ.copy()
        try:
            import vapoursynth

            pkgconfig_dir = Path(vapoursynth.__file__).parent / "pkgconfig"
            if pkgconfig_dir.is_dir():
                env["PKG_CONFIG_PATH"] = os.pathsep.join(
                    filter(None, [str(pkgconfig_dir), env.get("PKG_CONFIG_PATH", "")])
                )
                # pkgconf's --define-prefix (on by default on Windows) assumes <prefix>/lib/pkgconfig/*.pc
                # and breaks the wheel's vapoursynth/pkgconfig/*.pc include paths
                env["PKG_CONFIG_DONT_DEFINE_PREFIX"] = "1"
        except ImportError:
            pass

        setup_args = []
        # windows has nothing like auditwheel/delocate bundling GSL here so link it statically
        if sys.platform == "win32":
            setup_args.append("-Dstatic_gsl=true")

        subprocess.run([sys.executable, "-m", "mesonbuild.mesonmain", "setup", "build", *setup_args], check=True, env=env)
        subprocess.run([sys.executable, "-m", "mesonbuild.mesonmain", "compile", "-C", "build"], check=True)

        self.target_dir.mkdir(parents=True, exist_ok=True)
        for file_path in self.source_dir.glob("*"):
            if file_path.is_file() and file_path.suffix in [".dll", ".so", ".dylib"]:
                shutil.copy2(file_path, self.target_dir)

    def finalize(self, version: str, build_data: dict[str, Any], artifact_path: str) -> None:
        shutil.rmtree(self.target_dir.parent, ignore_errors=True)
