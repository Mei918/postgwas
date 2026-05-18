from __future__ import annotations

from dataclasses import dataclass, field
import os
from pathlib import Path
import shutil
import subprocess
from typing import Optional, Union


@dataclass
class ToolResult:
    command: list[str]
    returncode: int
    stdout: str
    stderr: str


@dataclass
class ExternalTool:
    executable: str
    extra_env: dict[str, str] = field(default_factory=dict)

    def resolve_executable(self) -> str:
        resolved = shutil.which(self.executable)
        if resolved is None:
            raise FileNotFoundError(
                f"Required executable not found on PATH: {self.executable}"
            )
        return resolved

    def run(self, args: list[str], cwd: Optional[Union[str, Path]] = None) -> ToolResult:
        command = [self.resolve_executable(), *args]
        env = os.environ.copy()
        env.update(self.extra_env)
        completed = subprocess.run(
            command,
            cwd=str(cwd) if cwd else None,
            capture_output=True,
            text=True,
            env=env,
            check=False,
        )
        return ToolResult(
            command=command,
            returncode=completed.returncode,
            stdout=completed.stdout,
            stderr=completed.stderr,
        )
