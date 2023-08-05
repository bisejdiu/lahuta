Lahuta models are simply classes which inherit from `Universe` and define fields as annotated attributes.

::: lahuta.core.universe.Universe
    options:
        show_root_heading: true
        merge_init_into_class: false
        group_by_category: false
        # explicit members list so we can set order and include `__init__` easily
        members:
          - __init__
          - ready
          - compute_neighbors
          - to

<!-- ::: pydantic.create_model
    options:
        show_root_heading: true -->