How to Edit Documentation
============================

The following is a basic guide for editing the documentation for this project. The documentation is rendered using Sphinx_ documentation builder

Quick start up
----------------------------------

#. Install Sphinx in Mac OS  ``brew install sphinx-doc`` or Linux (Ubuntu) ``apt-get install python3-sphinx``

#. Go to the docs folder ``cd docs``

#. Change any ``.rst`` file you need to modify

#. To visualise the results locally do ``make dirhtml`` and check the html files in the ``_build/dirhtml`` directory

#. When you are done, push your changes to the git repo.


Sphinx Workspace Structure
----------------------------

All of the files concerned with documentation can be found in the ``docs`` directory. At the top level there is the ``conf.py`` file and an ``index.rst`` file among others. There should be no need to change the ``conf.py`` file. The ``index.rst`` file is known as the master document and defines the document structure of the documentation (i.e. Menu Or Table of Contents structure). It contains the root of the “table of contents" tree -or toctree- that is used to connect the multiple files to a single hierarchy of documents. The TOC is defined using the ``toctree`` directive which is used as follows::

    .. toctree::
       :maxdepth: 2
       :caption: Getting Started

        usage/introduction
        usage/installation

The ``toctree`` inserts a TOC tree at the current location using the individual TOCs of the documents given in the directive command body. In other words if there are ``toctree`` directives in the files listed in the above example it will also be applied to the resulting TOC. Relative document names (not beginning with a slash) are relative to the document the directive occurs in, absolute names are relative to the source directory. Thus in the example above the ``usage`` directory is relative to the ``index.rst`` page . The ``:maxdepth:`` parameter defines the depth of the tree for that particular menu. The ``caption`` parameter is used to give a caption for that menu tree at that level. It should be noted the titles for the links of the menu items under that header would be taken from the titles of the referenced document. For example the menu item title for ``usage/introduction`` is taken from the main header specified in ``introduction.rst`` document in the ``usage`` directory. Also note the document name does not include the extention (i.e. .rst).

Thus the directory structure for the above example is shown below::

    ├── index.rst
    └── usage
        ├── introduction.rst
        └── installation.rst


Basic reStructuredText Syntax
-------------------------------

Now we will look at some basic reStructuredText syntax necessary to start editing the .rst files that are used to generate documentation. 

Headers
""""""""

**Section Header**

The following was used to make the header at the top of this page:
::

    How to Edit Documentation
    ==========================

**Subsection Header**

The follwoing was used to create the secondary header (e.g. Sphinx Workspace Structure section header)
::

    Sphinx Workspace structure
    ----------------------------

..... 


Lists
""""""
**Bullets List**
::

    - This is a bullet
    - This is a bullet

Will produce the following:

- This is a bullet
- This is a bullet


**Numbered List**
::

    #. This is a numbered list item
    #. This is a numbered list item

Will produce the following:

#. This is a numbered list item
#. This is a numbered list item

.....

Inline Markup
""""""""""""""
**Emphasis/Italics**
::

    *This is for emphasis*

Will produce the following 

*This is for emphasis*


**Bold**
::

    **This is bold text**

Will produce the following

**This is bold text**

..... 

**Code Sample**
::
    
    ``Backquotes = code sample``

Will produce the following:

``Backquotes = code sample``

**Apostraphies in Text**
::

    `don't know`

Will produce the following

`don't know`


**Literal blocks**

Literal code blocks are introduced by ending a paragraph with the special marker ``::``. The literal block must be indented (and, like all paragraphs, separated from the surrounding ones by blank lines)::

    This is a normal text paragraph. The next paragraph is a code sample::

        It is not processed in any way, except
        that the indentation is removed.

        It can span multiple lines.

    This is a normal text paragraph again.


The following is produced:

.....

This is a normal text paragraph. The next paragraph is a code sample::

    It is not processed in any way, except
    that the indentation is removed.

    It can span multiple lines.

This is a normal text paragraph again.

.....

**Doctest blocks**

Doctest blocks are interactive Python sessions cut-and-pasted into docstrings. They do not require the literal blocks syntax. The doctest block must end with a blank line and should not end with with an unused prompt:

>>> 1 + 1
2

**External links**

Use ```Link text <https://domain.invalid/>`_`` for inline web links `Link text <https://domain.invalid/>`_. If the link text should be the web address, you don’t need special markup at all, the parser finds links and mail addresses in ordinary text. *Important:* There must be a space between the link text and the opening ``<`` for the URL.

You can also separate the link and the target definition , like this
::

    This is a paragraph that contains `a link`_.

    .. _a link: https://domain.invalid/


Will produce the following:

This is a paragraph that contains `a link`_.

.. _a link: https://domain.invalid/



**Internal links**

Internal linking is done via a special reST role provided by Sphinx to cross-reference arbitrary locations. For this to work label names must be unique throughout the entire documentation. There are two ways in which you can refer to labels:

- If you place a label directly before a section title, you can reference to it with ``:ref:`label-name```. For example::

    .. _my-reference-label:

    Section to cross-reference
    --------------------------

    This is the text of the section.

    It refers to the section itself, see :ref:`my-reference-label`.

The ``:ref:`` role would then generate a link to the section, with the link title being “Section to cross-reference”. This works just as well when section and reference are in different source files. The above produces the following:

.....

.. _my-reference-label:

Section to cross-reference
"""""""""""""""""""""""""""

This is the text of the section.

It refers to the section itself, see :ref:`my-reference-label`.

.....

- Labels that aren’t placed before a section title can still be referenced, but you must give the link an explicit title, using this syntax: ``:ref:`Link title <label-name>```.


**Comments**

Every explicit markup block which isn’t a valid markup construct is regarded as a comment. For example::

    .. This is a comment.

Go to Sphinx_ for more documentation. 

.. _Sphinx: https://www.sphinx-doc.org
.. _reStructuredText: https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html

