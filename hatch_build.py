"""Build the web interface into src/scnascope/static when packaging, if it isn't built yet.

Needs Node.js (npm) at build time; the conda environment and recipe provide it. With a
prebuilt interface (src/scnascope/static/index.html present) nothing happens.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from hatchling.builders.hooks.plugin.interface import BuildHookInterface


class WebBuildHook(BuildHookInterface):
    PLUGIN_NAME = "web"

    def initialize(self, version, build_data):
        root = Path(self.root)
        static, web = root / "src" / "scnascope" / "static", root / "web"
        if (static / "index.html").exists() or not (web / "package.json").exists():
            return
        npm = shutil.which("npm")
        if npm is None:
            self.app.display_warning(
                "scnascope: Node.js (npm) not found, so the web interface was not built; "
                "the viewer will only show a notice. Install nodejs and reinstall."
            )
            return
        self.app.display_info("scnascope: building the web interface")
        if not (web / "node_modules").exists():
            subprocess.run([npm, "ci", "--no-audit", "--no-fund"], cwd=web, check=True)
        subprocess.run([npm, "run", "build"], cwd=web, check=True)
