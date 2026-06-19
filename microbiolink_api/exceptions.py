#!/usr/bin/env python

"""Custom exceptions for the user-facing MicrobioLink API."""


class MicrobioLinkAPIError(Exception):
    """Base exception for the new library-style MicrobioLink API."""


class InputFormatError(MicrobioLinkAPIError):
    """Raised when an input file or in-memory object has an invalid format."""
