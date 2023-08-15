Visualizing computed contacts is also easy and straighforward, as the following example shows:

???+ example "Example - Contact Map"
    ```py
    from lahuta import Luni

    luni = Luni("path/to/file.pdb")
    ns = luni.compute_neighbors()

    ns.plot() # (1)!
    ns.plot('full', half_only=True) # (2)!
    
    ```

    1. Plots a contact map for only the indices that are in the `pairs` array
    2. Plots a contact map for all atoms in the loaded system

## TODO:
Add example images