#!/usr/bin/env python

"""Custom exceptions for the user-facing MicrobioLink API."""


class MicrobioLinkError(Exception):
    """Base exception for the new library-style MicrobioLink API."""


class InputFormatError(MicrobioLinkError):
    """Raised when an input file or in-memory object has an invalid format."""
