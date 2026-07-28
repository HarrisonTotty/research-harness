import sys
import textwrap
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "docs"))

from lean_parser import (
    LeanDeclaration,
    LeanModule,
    LeanProse,
    dedent_block,
    parse_lean_source,
)

REPO_ROOT = Path(__file__).resolve().parents[1]


def parse(source: str) -> LeanModule:
    return parse_lean_source("Test", textwrap.dedent(source).lstrip("\n"))


def declarations(source: str) -> list[LeanDeclaration]:
    return [e for e in parse(source).entries if isinstance(e, LeanDeclaration)]


def only(source: str) -> LeanDeclaration:
    found = declarations(source)
    assert len(found) == 1
    return found[0]


# --------------------------------------------------------------------------- #
# Comment lexing
# --------------------------------------------------------------------------- #
def test_copyright_header_is_not_the_module_docstring():
    module = parse("""
        /-
        Copyright (c) 2026 Someone. All rights reserved.
        -/
        import Mathlib.Tactic

        /-!
        # Real module docstring
        -/
    """)
    assert module.docstring == "# Real module docstring"


def test_module_without_docstring():
    assert parse("def f : Nat := 0").docstring is None


def test_nested_block_comments_do_not_close_early():
    source = """
        /-!
        # Title

        Mentions a /- nested -/ comment.
        -/
        def f : Nat := 0
    """
    module = parse(source)
    assert module.docstring is not None
    assert "nested" in module.docstring
    # The inner `-/` must not end the docstring and leave its tail as code.
    assert [d.name for d in declarations(source)] == ["f"]


def test_later_module_blocks_become_prose_in_source_order():
    module = parse("""
        /-! # Header -/

        def a : Nat := 0

        /-! ## Section two -/

        def b : Nat := 1
    """)
    kinds = [type(entry).__name__ for entry in module.entries]
    assert kinds == ["LeanDeclaration", "LeanProse", "LeanDeclaration"]
    prose = module.entries[1]
    assert isinstance(prose, LeanProse)
    assert prose.body == "## Section two"


def test_docstring_markers_are_stripped_from_bodies():
    declaration = only("""
        /-- The doc. -/
        def f : Nat := 0
    """)
    assert declaration.docstring == "The doc."


def test_multiline_docstring_is_dedented():
    declaration = only("""
        /-- First line.

        Second paragraph, indented in source.
        -/
        def f : Nat := 0
    """)
    assert declaration.docstring == (
        "First line.\n\nSecond paragraph, indented in source."
    )


# --------------------------------------------------------------------------- #
# Declarations
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    ("source", "kind", "name", "signature"),
    [
        (
            "def double (n : Nat) : Nat := n + n",
            "def",
            "double",
            "def double (n : Nat) : Nat",
        ),
        (
            "theorem t (n : Nat) : n = n := by rfl",
            "theorem",
            "t",
            "theorem t (n : Nat) : n = n",
        ),
        ("lemma l : True := trivial", "lemma", "l", "lemma l : True"),
        ("abbrev A := Nat", "abbrev", "A", "abbrev A"),
        ("axiom ax : True", "axiom", "ax", "axiom ax : True"),
        ("class C (a : Type) where", "class", "C", "class C (a : Type)"),
        ("structure S where", "structure", "S", "structure S"),
        ("inductive I where", "inductive", "I", "inductive I"),
    ],
)
def test_declaration_kinds_and_signatures(source, kind, name, signature):
    declaration = only(source)
    assert (declaration.kind, declaration.name, declaration.signature) == (
        kind,
        name,
        signature,
    )


def test_modifiers_are_kept_in_the_signature_but_not_the_name():
    declaration = only("noncomputable def f : Nat := 0")
    assert declaration.name == "f"
    assert declaration.signature == "noncomputable def f : Nat"


def test_signature_spanning_multiple_lines_is_joined():
    declaration = only("""
        theorem long (a : Nat)
            (b : Nat) :
            a + b = b + a := by
          omega
    """)
    assert declaration.signature == "theorem long (a : Nat) (b : Nat) : a + b = b + a"


def test_where_ends_a_signature():
    declaration = only("""
        structure Point where
          x : Nat
          y : Nat
    """)
    assert declaration.signature == "structure Point"


def test_inductive_constructors_do_not_leak_into_the_signature():
    declaration = only("""
        inductive Color : Type
          | red
          | green
    """)
    assert declaration.signature == "inductive Color : Type"


def test_anonymous_declarations_are_skipped():
    assert declarations("instance : Inhabited Nat := ⟨0⟩") == []


def test_examples_are_skipped():
    assert declarations("example : True := trivial") == []


def test_undocumented_declarations_are_still_reported():
    declaration = only("def f : Nat := 0")
    assert declaration.docstring is None


def test_attributes_between_docstring_and_declaration_are_transparent():
    declaration = only("""
        /-- Documented. -/
        @[simp]
        def f : Nat := 0
    """)
    assert declaration.docstring == "Documented."


def test_docstring_does_not_leak_past_intervening_code():
    declaration = only("""
        /-- Belongs to nothing. -/
        variable (n : Nat)

        def f : Nat := 0
    """)
    assert declaration.docstring is None


# --------------------------------------------------------------------------- #
# Namespaces
# --------------------------------------------------------------------------- #
def test_namespace_qualifies_declaration_names():
    declaration = only("""
        namespace Foo
        def f : Nat := 0
        end Foo
    """)
    assert declaration.name == "Foo.f"


def test_nested_namespaces_compose():
    declaration = only("""
        namespace Foo
        namespace Bar
        def f : Nat := 0
        end Bar
        end Foo
    """)
    assert declaration.name == "Foo.Bar.f"


def test_dotted_namespace_is_kept_whole():
    declaration = only("""
        namespace Foo.Bar
        def f : Nat := 0
        end Foo.Bar
    """)
    assert declaration.name == "Foo.Bar.f"


def test_namespace_closes_so_later_declarations_are_unqualified():
    found = declarations("""
        namespace Foo
        def inside : Nat := 0
        end Foo
        def outside : Nat := 0
    """)
    assert [d.name for d in found] == ["Foo.inside", "outside"]


def test_section_end_does_not_pop_the_enclosing_namespace():
    declaration = only("""
        namespace Foo
        section Helpers
        end Helpers
        def f : Nat := 0
        end Foo
    """)
    assert declaration.name == "Foo.f"


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def test_dedent_block_on_empty_body():
    assert dedent_block("\n   \n") == ""


def test_dedent_block_keeps_relative_indentation():
    """A body starting on its own line dedents by its true common prefix."""
    body = "\n    Intro.\n\n        code block\n    "
    assert dedent_block(body) == "Intro.\n\n    code block"


def test_dedent_block_when_first_line_shares_the_delimiter():
    """The first line has no indent of its own, so it must not set the prefix."""
    body = " Intro.\n\n    Body.\n\n        code block\n"
    assert dedent_block(body) == "Intro.\n\nBody.\n\n    code block"


# --------------------------------------------------------------------------- #
# The repository's own Lean sources
# --------------------------------------------------------------------------- #
def test_parses_the_repository_lean_library():
    sources = sorted((REPO_ROOT / "src" / "theorems").rglob("*.lean"))
    assert sources, "expected at least one Lean source to document"
    for source in sources:
        module = parse_lean_source(source.stem, source.read_text(encoding="utf-8"))
        assert module.docstring, f"{source} has no module docstring"


def test_basic_module_declarations_round_trip():
    source = REPO_ROOT / "src" / "theorems" / "Theorems" / "Basic.lean"
    module = parse_lean_source("Theorems.Basic", source.read_text(encoding="utf-8"))
    found = {
        entry.name: entry
        for entry in module.entries
        if isinstance(entry, LeanDeclaration)
    }
    assert found["ResearchHarness.double"].signature == "def double (n : ℕ) : ℕ"
    assert found["ResearchHarness.double_eq_two_mul"].signature == (
        "theorem double_eq_two_mul (n : ℕ) : double n = 2 * n"
    )
