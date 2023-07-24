import numpy as np
from numpy.typing import NDArray


def calc_pairwise_distances(matrix1: NDArray[np.float32], matrix2: NDArray[np.float32]) -> NDArray[np.float32]:
    """
    Calculate the pairwise Euclidean distance between two sets of vectors.

    This function computes the Euclidean distance between each pair of vectors,
    where one vector is taken from `matrix1` and each 3D vector from the second dimension of `matrix2`.
    The output is a 2D array where the entry at position (i, j) is the distance between the i-th vector
    in `matrix1` and the j-th vector in `matrix2`.

    Args:
        matrix1 (NDArray[np.float32]): Shape (n, 3). The first set of vectors. Each row represents a vector.
        matrix2 (NDArray[np.float32]): Shape (n, m, 3). The second set of vectors.
                                                        Each row in the second dimension represents a vector.

    Returns:
        NDArray[np.float32]: Shape (n, m). The pairwise distance matrix.


    Example:
        >>> mat1 = np.array([[1, 0, 0], [0, 1, 0]])
        >>> mat2 = np.array([[[0, 0, 0], [1, 1, 1]], [[1, 1, 1], [0, 0, 0]]])
        >>> calc_pairwise_distances(mat1, mat2)
        array([[[1.        , 1.41421356],
                [1.41421356, 1.        ]]])

    """
    reshaped_matrix1 = matrix1[:, np.newaxis, :]  # type: ignore
    distance_matrix = np.linalg.norm(reshaped_matrix1 - matrix2, axis=-1)

    return distance_matrix


# pylint: disable=W1114
def calc_vertex_angles(
    vertex: NDArray[np.float32],
    point1: NDArray[np.float32],
    point2: NDArray[np.float32],
    degrees: bool = False,
) -> NDArray[np.float32]:
    """Calculate the angle between three sets of points.

    This function calculates the angle at each `vertex` created by points `point1` and `point2`.

    Args:
        vertex (NDArray[np.float32]): Shape (n, m, 3). The coordinates of the vertex point(s), where the angle is being measured.
        point1 (NDArray[np.float32]): Shape (n, m, 3). The coordinates of the first point(s).
        point2 (NDArray[np.float32]): Shape (n, m, 3). The coordinates of the second point(s).
        degrees (bool, optional): If True, the angle is returned in degrees. Otherwise, it's returned in radians.
            Default is False (radians).

    Returns:
        NDArray[np.float32]: Shape (n, m). The calculated angle(s) in radians (default) or degrees,
                                            depending on the value of the `degrees` argument.

    """
    # Vector from vertex to point1 & point2
    vector1: NDArray[np.float32] = point1 - vertex
    vector2: NDArray[np.float32] = point2 - vertex

    # Normalize vector1
    magnitude_vector1: NDArray[np.float32] = np.linalg.norm(vector1, axis=-1)
    normalized_vector1: NDArray[np.float32] = vector1 / magnitude_vector1[:, :, np.newaxis]

    # Normalize vector2
    magnitude_vector2: NDArray[np.float32] = np.linalg.norm(vector2, axis=-1)
    normalized_vector2: NDArray[np.float32] = vector2 / magnitude_vector2[:, :, np.newaxis]

    # Dot product of normalized vectors
    dot_product: NDArray[np.float32] = np.sum(normalized_vector1 * normalized_vector2, axis=-1)

    angle_rad: NDArray[np.float32] = np.arccos(dot_product)

    if degrees:
        return np.degrees(angle_rad)

    return angle_rad


def calc_vec_angle(vector1: NDArray[np.float32], vector2: NDArray[np.float32]) -> NDArray[np.float32]:
    """
    Calculates the angle between two vectors in degrees.

    This function computes the angle between pairs of vectors. Each vector in the pair is normalized first,
    and then the dot product is calculated. The angle is obtained using the arccosine of the dot product,
    adjusted by the sign of the dot product. Finally, the angle is converted to degrees.

    Args:
        vector1 (NDArray[np.float32]): The first vector or a set of vectors in each pair.
        vector2 (NDArray[np.float32]): The second vector or a set of vectors in each pair.

    Returns:
        NDArray[np.float32]: The angles between pairs of vectors in degrees.

    Example:
        >>> v1 = np.array([[1, 0, 0], [0, 1, 0]])
        >>> v2 = np.array([[0, 1, 0], [0, 0, 1]])
        >>> calc_vec_angle(v1, v2)
        array([90., 90.], dtype=float32)
    """
    normalized_vector1 = vector1 / np.linalg.norm(vector1, axis=-1, keepdims=True)
    normalized_vector2 = vector2 / np.linalg.norm(vector2, axis=-1, keepdims=True)
    dot_product = np.einsum("ij,ij->i", normalized_vector1, normalized_vector2)  # type: ignore
    raw_angle = np.arccos(np.clip(dot_product, -1.0, 1.0))
    adjusted_angle = np.sign(dot_product) * raw_angle
    angle_in_degrees = np.degrees(adjusted_angle)
    angle_in_degrees[angle_in_degrees < 0] = 180 + angle_in_degrees[angle_in_degrees < 0]

    return angle_in_degrees


def calc_vec_line_angles(vector: NDArray[np.float_], line_direction: NDArray[np.float_]) -> NDArray[np.float_]:
    """
    Calculate the angle between a vector and a line direction in degrees.

    This function computes the angle between each pair of vector and line direction.
    The dot product of the vector and the normalized line direction is calculated.
    The raw angle is obtained using the arccosine of the dot product, and it is adjusted by the sign of the dot product.
    If the adjusted angle is less than zero, it is incremented by pi.
    Finally, the angle is converted to degrees.

    Args:
        vector (NDArray[np.float_]): The vector or set of vectors.
        line_direction (NDArray[np.float_]): The direction or set of directions of the lines.

    Returns:
        NDArray[np.float_]: The angles between the vectors and line directions in degrees.

    Example:
        >>> vec = np.array([[1, 0, 0], [0, 1, 0]])
        >>> line_dir = np.array([[0, 1, 0], [0, 0, 1]])
        >>> calc_vec_line_angles(vec, line_dir)
        array([90., 90.], dtype=float32)
    """
    normalized_line_direction = line_direction / np.linalg.norm(line_direction, axis=1)[..., np.newaxis]
    dot_product = np.einsum("ij,ij->i", vector, normalized_line_direction)  # type: ignore
    raw_angle = np.arccos(dot_product)
    adjusted_angle = np.sign(dot_product) * raw_angle

    adjusted_angle = np.where(adjusted_angle < 0, adjusted_angle + np.pi, adjusted_angle)

    return np.degrees(adjusted_angle)
