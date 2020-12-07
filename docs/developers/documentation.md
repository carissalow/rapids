# Documentation

We use [mkdocs](https://www.mkdocs.org/) with the [material theme](https://squidfunk.github.io/mkdocs-material/) to write these docs. Whenever you make any changes, just push them back to the repo and the documentation will be deployed automatically.

## Set up development environment

1. Make sure your conda environment is active
2. `pip install mkdocs`
3. `pip install mkdocs-material`

## Preview

Run the following command in RAPIDS root folder and go to [http://127.0.0.1:8000](http://127.0.0.1:8000):

```bash
mkdocs serve
```

## File Structure

The documentation config file is `/mkdocs.yml`, if you are adding new `.md` files to the docs modify the `nav` attribute at the bottom of that file. You can use the hierarchy there to find all the files that appear in the documentation.

## Reference

Check this [page](https://squidfunk.github.io/mkdocs-material/reference/abbreviations/) to get familiar with the different visual elements we can use in the docs (admonitions, code blocks, tables, etc.) You can also refer to `/docs/setup/installation.md` and `/docs/setup/configuration.md` to see practical examples of these elements.

!!! hint
    Any links to internal pages should be relative to the current page. For example, any link from this page (documentation) which is inside `./developers` should begin with `../` to go one folder level up like:
    ```md
    [mylink](../setup/installation.md)
    ```

## Extras

You can insert [emojis](https://facelessuser.github.io/pymdown-extensions/extensions/emoji/) using this syntax `:[SOURCE]-[ICON_NAME]` from the following sources:

- https://materialdesignicons.com/
- https://fontawesome.com/icons/tasks?style=solid
- https://primer.style/octicons/

You can use this [page](https://www.tablesgenerator.com/markdown_tables) to create markdown tables more easily
