"""Build hooks for the documentation site.

Registered through the ``hooks`` key in ``properdocs.yml``, which appends them
after every plugin, so these run once the plugins they depend on have had their
turn.

Imports come from ``properdocs``, the generator that actually runs the build,
rather than from ``mkdocs``, which several plugins still pull in. The two ship
parallel copies of the same classes: ``InclusionLevel`` in particular is a
distinct enum in each, and members of one do not compare equal to members of the
other, so a value taken from the wrong package would silently fail the ``==``
test ProperDocs applies to it.
"""

from collections.abc import Mapping
from typing import cast

from properdocs.config.defaults import ProperDocsConfig
from properdocs.structure.files import Files, InclusionLevel
from properdocs.structure.nav import Navigation

_LITERATE_NAV = "literate-nav"
_DEFAULT_NAV_FILE = "SUMMARY.md"


def _nav_file_name(config: ProperDocsConfig) -> str:
    """Return the file name ``literate-nav`` is configured to read navs from.

    Read from the plugin rather than assumed, so changing ``nav_file`` in
    ``properdocs.yml`` does not silently leave this hook matching nothing.
    """
    # `PluginCollection` is loosely typed upstream; only the name is needed.
    plugins = cast("Mapping[str, object]", config.plugins)
    plugin_config: object = getattr(plugins.get(_LITERATE_NAV), "config", None)
    nav_file: object = getattr(plugin_config, "nav_file", None)
    return nav_file if isinstance(nav_file, str) else _DEFAULT_NAV_FILE


def on_nav(nav: Navigation, /, *, config: ProperDocsConfig, files: Files) -> Navigation:
    """Drop ``literate-nav``'s nav files from the built site.

    The plugin marks them ``NOT_IN_NAV``, which keeps them out of the menu but
    still renders each one as a page — leaving stub pages in the output, the
    sitemap, and the search index, and emitting link warnings for the directory
    links they are made of. They are scaffolding, not documentation, so they are
    excluded outright.

    This runs on ``on_nav`` rather than ``on_files`` because ``literate-nav``
    needs to read the files first, and pages to render are selected after the
    nav event.
    """
    nav_file = _nav_file_name(config)

    for file in files:
        if file.src_uri == nav_file or file.src_uri.endswith(f"/{nav_file}"):
            file.inclusion = InclusionLevel.EXCLUDED

    return nav
