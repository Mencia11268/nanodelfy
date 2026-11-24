"""El servidor web debe poder importarse con las dependencias declaradas.

Regresión: jinja2 faltaba en pyproject.toml y `nanodelfy web` fallaba en
imágenes Docker limpias (ImportError de Jinja2Templates), aunque en
desarrollo local funcionaba porque jinja2 ya estaba presente.
"""
import subprocess
import sys


def test_server_module_imports():
    from nanodelfy.server import app  # noqa: F401
    assert app.title


def test_web_entrypoint_importable_in_clean_env():
    """Simula el arranque real: importa el módulo CLI que usa `nanodelfy web`."""
    code = "from nanodelfy.cli import _cmd_web; from nanodelfy.server import app"
    subprocess.run([sys.executable, "-c", code], check=True)
