+++
title = "Guidelines for contributions to the documentation"
author = "Feliks Kiszkurno"
weight = 1023
+++

## Documentation guidelines

1. Providing the documentation for new and advanced users is part of the development.
2. Keep the entry level low - require as little preexisting knowledge from users as possible
3. Provide context - Link to other parts of documentation.
   It will provide context for new users and help them discover features they were not aware of.
4. Be explicit - Things that seem self explanatory to you, may not be obvious to others.
5. Use templates - This will help readers find what they need faster and makes writing docs easier.
6. Completeness - describe all features of what you implemented.
   The code is not sufficient documentation on its own.

## Technical recommendations

1. In the Markdown files, each sentence should be placed in a separate line.
It will make reviewing and tracking of changes easier.
2. Missing content can be marked by adding `TODO:` followed by a brief description of the missing content.
The line containing `TODO:` should be placed below the heading of the section it refers to.

## Style recommendations

1. Linking directly to the relevant paragraph should be preferred over linking to the whole page (see code snippet below point 2.).
2. Titles of links to other pages should be written in bold:
```md
[**Title of the link**](/link/to/page/#paragraph)
```

## Links to templates

[Process description](https://gitlab.opengeosys.org/ogs/ogs/-/tree/master/web/content/docs/devguide/documentation/docs-guidelines/templates/process/index.md)
