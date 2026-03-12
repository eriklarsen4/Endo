# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 09:37:10 2026

@author: Erik
"""

# %% Imports

import importlib
import os
import pkgutil
import saint


# %% Helpers

def _module_exists(module_path):
    try:
        importlib.import_module(module_path)
        return True
    except Exception:
        return False


# %% Tests

def test_all_expected_modules_importable():
    """
    Ensure that all core modules in the SAINT package import cleanly.
    This catches missing files, missing __init__.py files, and broken imports.
    """

    expected_modules = [
        "saint.io.data_input",
        "saint.pipeline.classical_saint",
        "saint.pipeline.hierarchical_saint",
        "saint.model.classical_em_wrapper",
        "saint.model.hierarchical_em_wrapper",
        "saint.model.classical_responsibilities",
        "saint.model.hierarchical_responsibilities",
        "saint.model.updates.lambda_updates",
        "saint.model.updates.pi_updates",
        "saint.model.updates.tau_updates",
        "saint.model.classical_likelihood",
        "saint.model.hierarchical_likelihood",
        "saint.diagnostics.diagnostics_classical",
        "saint.diagnostics.diagnostics_hierarchical",
    ]

    for module in expected_modules:
        assert _module_exists(module), f"Module failed to import: {module}"


def test_all_directories_have_init_files():
    """
    Ensure every directory in the saint package contains an __init__.py file.
    This enforces stable imports.
    """

    saint_path = os.path.dirname(saint.__file__)

    for root, dirs, files in os.walk(saint_path):
        # Skip hidden directories like __pycache__
        dirs[:] = [d for d in dirs if not d.startswith("__")]

        if root.endswith("saint"):
            # top-level package already guaranteed to have __init__.py
            continue

        assert "__init__.py" in files, f"Missing __init__.py in {root}"


def test_no_broken_submodules():
    """
    Ensure that pkgutil can walk the entire saint package without errors.
    This catches namespace issues and missing modules.
    """

    for module_info in pkgutil.walk_packages(saint.__path__, prefix="saint."):
        module_name = module_info.name
        assert _module_exists(module_name), f"Broken submodule: {module_name}"