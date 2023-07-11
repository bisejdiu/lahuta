def get_class_methods(cls):
    """Returns a list of class methods."""
    return [attr for attr in dir(cls) if callable(getattr(cls, attr))]


def get_class_properties(cls):
    """Returns a list of class properties."""
    return [attr for attr in dir(cls) if isinstance(getattr(cls, attr), property)]


def get_class_attributes(cls):
    """Returns a list of class attributes that are not callable."""
    # properties = get_class_properties(cls)
    # print("properties", properties)
    return [
        attr
        for attr in dir(cls)
        if not attr.startswith("__") and not callable(getattr(cls, attr))
        # and attr not in properties
    ]
