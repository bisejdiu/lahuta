from __future__ import annotations as _annotations

import json
import logging
import os
import re
from pathlib import Path
from textwrap import indent

from mkdocs.config import Config
from mkdocs.structure.files import Files
from mkdocs.structure.pages import Page

logger = logging.getLogger('mkdocs.plugin')
THIS_DIR = Path(__file__).parent
DOCS_DIR = THIS_DIR.parent
PROJECT_ROOT = DOCS_DIR.parent


def on_pre_build(config: Config) -> None:
    """
    Before the build starts.
    """
    add_changelog()


def on_files(files: Files, config: Config) -> Files:
    """
    After the files are loaded, but before they are read.
    """
    return files


def on_page_markdown(markdown: str, page: Page, config: Config, files: Files) -> str:
    """
    Called on each file after it is read and before it is converted to HTML.
    """
    markdown = insert_json_output(markdown)
    markdown = remove_code_fence_attributes(markdown)
    if md := render_index(markdown, page):
        return md
    else:
        return markdown


def add_changelog() -> None:
    history = (PROJECT_ROOT / 'HISTORY.md').read_text()
    history = re.sub(r'(\s)@([\w\-]+)', r'\1[@\2](https://github.com/\2)', history, flags=re.I)
    history = re.sub(r'\[GitHub release]\(', r'[:simple-github: GitHub release](', history)
    history = re.sub('@@', '@', history)
    new_file = DOCS_DIR / 'changelog.md'

    # avoid writing file unless the content has changed to avoid infinite build loop
    if not new_file.is_file() or new_file.read_text() != history:
        new_file.write_text(history)


MIN_MINOR_VERSION = 10
MAX_MINOR_VERSION = 12


def insert_json_output(markdown: str) -> str:
    """
    Find `output="json"` code fence tags and replace with a separate JSON section
    """

    def replace_json(m: re.Match[str]) -> str:
        start, attrs, code = m.groups()

        def replace_last_print(m2: re.Match[str]) -> str:
            ind, json_text = m2.groups()
            json_text = indent(json.dumps(json.loads(json_text), indent=2), ind)
            # no trailing fence as that's not part of code
            return f'\n{ind}```\n\n{ind}JSON output:\n\n{ind}```json\n{json_text}\n'

        code = re.sub(r'\n( *)"""(.*?)\1"""\n$', replace_last_print, code, flags=re.S)
        return f'{start}{attrs}{code}{start}\n'

    return re.sub(r'(^ *```)([^\n]*?output="json"[^\n]*?\n)(.+?)\1', replace_json, markdown, flags=re.M | re.S)


def remove_code_fence_attributes(markdown: str) -> str:
    """
    There's no way to add attributes to code fences that works with both pycharm and mkdocs, hence we use
    `py key="value"` to provide attributes to pytest-examples, then remove those attributes here.

    https://youtrack.jetbrains.com/issue/IDEA-297873 & https://python-markdown.github.io/extensions/fenced_code_blocks/
    """

    def remove_attrs(match: re.Match[str]) -> str:
        suffix = re.sub(
            r' (?:test|lint|upgrade|group|requires|output|rewrite_assert)=".+?"', '', match.group(2), flags=re.M
        )
        return f'{match.group(1)}{suffix}'

    return re.sub(r'^( *``` *py)(.*)', remove_attrs, markdown, flags=re.M)


def render_index(markdown: str, page: Page) -> str | None:
    if page.file.src_uri != 'index.md':
        return None

    if version := os.getenv('LAHUTA_VERSION'):
        url = f'https://github.com/bisejdiu/lahuta/releases/tag/{version}'
        version_str = f'Documentation for version: [{version}]({url})'
    elif (version_ref := os.getenv('GITHUB_REF')) and version_ref.startswith('refs/tags/'):
        version = re.sub('^refs/tags/', '', version_ref.lower())
        url = f'https://github.com/bisejdiu/lahuta/releases/tag/{version}'
        version_str = f'Documentation for version: [{version}]({url})'
    elif sha := os.getenv('GITHUB_SHA'):
        url = f'https://github.com/bisejdiu/bisejdiu/commit/{sha}'
        sha = sha[:7]
        version_str = f'Documentation for development version: [{sha}]({url})'
    else:
        version_str = 'Documentation for development version'
    logger.info('Setting version prefix: %r', version_str)
    return re.sub(r'{{ *version *}}', version_str, markdown)
