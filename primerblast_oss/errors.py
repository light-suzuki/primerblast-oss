"""Structured exception hierarchy for the public library API.

Every class subclasses :class:`RuntimeError` so legacy code that catches
``RuntimeError`` keeps working while new code can branch on the specific
type. See ``docs/PUBLIC_API.md`` for the compatibility policy.
"""
from __future__ import annotations


class PrimerblastError(RuntimeError):
    """Base class for every primerblast-oss library error."""


class ToolMissingError(PrimerblastError):
    """A required external executable (primer3_core / blastn / makeblastdb)
    could not be found."""


class InvalidDatabaseError(PrimerblastError):
    """A BLAST database or genome index is missing or unusable."""


class Primer3Error(PrimerblastError):
    """The primer3_core process failed."""


class BlastError(PrimerblastError):
    """A blastn / makeblastdb process failed."""


class MalformedInputError(PrimerblastError):
    """Input sequences, coordinates, or options are malformed."""


class SearchIncompleteError(PrimerblastError):
    """Search evidence is incomplete (truncated or unparseable), so a
    definitive specificity statement is not possible."""


class CancelledError(PrimerblastError):
    """A long-running library call was cancelled by the caller's callback."""
