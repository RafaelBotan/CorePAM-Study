"""Inject DOIs from _REF_VERIFY.md into tese_references.bib.

Parses each `### N. bibkey — STATUS` section and its DOI line,
adds `doi = {...}` field to the corresponding bib entry if missing.
"""
import re

BIB = 'tese_references.bib'
VERIFY = '_REF_VERIFY.md'

with open(VERIFY, encoding='utf-8') as f:
    verify_text = f.read()

# Extract mapping: bibkey -> doi
key_to_doi = {}
for m in re.finditer(r'###\s*\d+\.\s*(\S+)(?:\s*\([^)]*\))?\s*—\s*.*?\n((?:(?!###\s).)*)', verify_text, re.S):
    key = m.group(1).strip()
    body = m.group(2)
    # Find first DOI line (with or without colon)
    doi_m = re.search(r'\bDOI[:\s]+([0-9][^\s—]+)', body)
    if doi_m:
        doi = doi_m.group(1).rstrip('.,;')
        key_to_doi[key] = doi

print(f'Extracted {len(key_to_doi)} DOIs from {VERIFY}')

with open(BIB, encoding='utf-8') as f:
    bib = f.read()

# Insert doi field in each entry that lacks one. Entry = @type{key, ... }
entries = list(re.finditer(r'@\w+\{([^,]+),\s*\n(.*?)\n\}', bib, re.S))
print(f'Found {len(entries)} bib entries')

new_bib_parts = []
last_end = 0
added = 0
skipped = 0
no_doi_found = []

for m in entries:
    key = m.group(1).strip()
    body = m.group(2)
    start, end = m.span()
    new_bib_parts.append(bib[last_end:start])
    full_entry = m.group(0)
    if re.search(r'\n\s*doi\s*=', full_entry, re.I):
        new_bib_parts.append(full_entry)
        skipped += 1
    else:
        doi = key_to_doi.get(key)
        if doi:
            # Insert "  doi = {...}" before the closing }
            # Strip trailing whitespace/newlines from body, add comma to last field, add doi line
            trimmed = full_entry.rstrip()
            if trimmed.endswith('}'):
                inner = trimmed[:-1].rstrip()
                # ensure inner last line ends in comma
                if not inner.endswith(','):
                    inner = inner + ','
                new_full = inner + f'\n  doi = {{{doi}}}\n}}'
                new_bib_parts.append(new_full)
                added += 1
            else:
                new_bib_parts.append(full_entry)
                skipped += 1
        else:
            new_bib_parts.append(full_entry)
            no_doi_found.append(key)
    last_end = end

new_bib_parts.append(bib[last_end:])
new_bib = ''.join(new_bib_parts)

with open(BIB, 'w', encoding='utf-8') as f:
    f.write(new_bib)

print(f'Added DOI to {added} entries')
print(f'Skipped (already had DOI) {skipped} entries')
print(f'No DOI found in _REF_VERIFY for: {no_doi_found}')
