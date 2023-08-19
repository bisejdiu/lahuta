"""Defines several mathematical functions for scientific computations, particularly
geometric calculations for molecular structures. The functions are built on the numpy 
library for efficient vectorized calculations.

Functions:
    ```
    - calc_pairwise_distances(matrix1, matrix2): Calculates pairwise Euclidean distance between two sets of vectors.
    - calc_vertex_angles(vertex, point1, point2, degrees=False): Calculates the angle between three sets of points.
    - calc_vec_angle(vector1, vector2): Calculates the angle between two vectors in degrees.
    - calc_vec_line_angles(vector, line_direction): \
        Calculates the angle between a vector and a line direction in degrees.
    - distance(vector1, vector2): Calculates the Euclidean distance between two vectors or sets of vectors.
    - normalize(vector): Normalizes a vector or a set of vectors.
    - dot_product(vector1, vector2): Calculates the dot product of two 1D vectors or sets of 1D vectors.
    - dot_product_3d(vector1, vector2): Calculates the dot product along the last axis of \
        two 2D vectors or sets of 2D vectors.
    ```
    
??? example "Example"
    ```py
    import numpy as np
    from math import calc_pairwise_distances, calc_vertex_angles, calc_vec_angle, calc_vec_line_angles

    v1 = np.array([[1, 0, 0], [0, 1, 0]])
    v2 = np.array([[0, 0, 0], [1, 1, 1]])
    v3 = np.array([[0, 0, 1], [1, 1, 0]])

    distances = calc_pairwise_distances(v1, v2)
    vertex_angles = calc_vertex_angles(v1, v2, v3, degrees=True)
    vec_angles = calc_vec_angle(v1, v2)
    vec_line_angles = calc_vec_line_angles(v1, v2)

    print(f"Distances: {distances}")
    print(f"Vertex Angles: {vertex_angles}")
    print(f"Vector Angles: {vec_angles}")
    print(f"Vector Line Angles: {vec_line_angles}")
    ```
"""

import numpy as np
from numpy.typing import NDArray


def distance(vector1: NDArray[np.float32], vector2: NDArray[np.float32]) -> NDArray[np.float32]:
    """Compute the Euclidean distance between two vectors or sets of vectors.

    This function computes the Euclidean distance between each pair of vectors,
    where each pair consists of one vector from `vector1` and one vector from `vector2`.
    The computation is performed along the last axis (i.e., `axis=-1`).

    Args:
        vector1 (NDArray[np.float32]): The first vector or set of vectors. Shape can be (n,) for a single vector,
            or (m, n) for a set of m vectors, or (k, m, n) for a set of k groups of m vectors.
        vector2 (NDArray[np.float32]): The second vector or set of vectors. The shape should match `vector1`.

    Returns:
        NDArray[np.float32]: Euclidean distance between each pair of vectors. The shape is (m,) for a set of m vectors,
            or (k, m) for a set of k groups of m vectors.

    ??? example "Example"
        ```py
        vec1 = np.array([[1, 0, 0], [0, 1, 0]])
        vec2 = np.array([[0, 1, 0], [0, 0, 1]])
        distance(vec1, vec2)
        array([1.41421356, 1.        ], dtype=float32)
        ```
    """
    result: NDArray[np.float32] = np.linalg.norm(vector1 - vector2, axis=-1)
    return result


def normalize(vector: NDArray[np.float32]) -> NDArray[np.float32]:
    """Normalize a vector or a set of vectors.

    This function normalizes a vector (or set of vectors) by dividing each vector
    by its corresponding L2 norm (Euclidean norm). The L2 norm is computed along the last axis.

    Args:
        vector (NDArray[np.float32]): Input vector or set of vectors. The shape can be (n,) for a single vector,
            or (m, n) for a set of m vectors, or (k, m, n) for a set of k groups of m vectors.

    Returns:
        NDArray[np.float32]: Normalized vector or set of vectors. The shape is the same as the input vector.

    ??? example "Example"
        ```py
        vec = np.array([[1, 0, 0], [0, 1, 0]])
        normalize(vec)
        array([[1., 0., 0.],
               [0., 1., 0.]], dtype=float32)
        ```
    """
    result: NDArray[np.float32] = vector / np.linalg.norm(vector, axis=-1, keepdims=True)
    return result


def dot_product(vector1: NDArray[np.float32], vector2: NDArray[np.float32]) -> NDArray[np.float32]:
    """Compute the dot product of two 1D vectors or sets of 1D vectors.

    This function uses Einstein summation notation to calculate the dot product.
    The result is a 1D array with the dot product of the corresponding pairs of 1D vectors.

    Args:
        vector1 (NDArray[np.float32]): First vector or set of vectors. Shape is (n,) or (m, n).
        vector2 (NDArray[np.float32]): Second vector or set of vectors. Shape is (n,) or (m, n).

    Returns:
        NDArray[np.float32]: The dot product of vector1 and vector2. Shape is () or (m,).
    """
    return np.einsum("ij,ij->i", vector1, vector2)  # type: ignore


def dot_product_3d(vector1: NDArray[np.float32], vector2: NDArray[np.float32]) -> NDArray[np.float32]:
    """Compute the dot product along the last axis of two 2D vectors or sets of 2D vectors.

    This function uses Einstein summation notation to calculate the dot product.
    The result is a 2D array with the dot product of the corresponding pairs of 2D vectors.

    Args:
        vector1 (NDArray[np.float32]): First vector or set of vectors. The shape is either (n, m) or (k, n, m).
        vector2 (NDArray[np.float32]): Second vector or set of vectors. The shape is either (n, m) or (k, n, m).

    Returns:
        NDArray[np.float32]: The dot product of vector1 and vector2 along the last axis. \
            The shape is either (n,) or (k, n).
    """
    return np.einsum("ijk,ijk->ij", vector1, vector2)  # type: ignore


def calc_pairwise_distances(matrix1: NDArray[np.float32], matrix2: NDArray[np.float32]) -> NDArray[np.float32]:
    """Calculate the pairwise Euclidean distance between two sets of vectors.

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


    ??? example "Example"
        ```py
        mat1 = np.array([[1, 0, 0], [0, 1, 0]])
        mat2 = np.array([[[0, 0, 0], [1, 1, 1]], [[1, 1, 1], [0, 0, 0]]])
        calc_pairwise_distances(mat1, mat2)
        array([[[1.        , 1.41421356],
                [1.41421356, 1.        ]]])
        ```
    """
    reshaped_matrix1 = matrix1[:, np.newaxis, :]  # type: ignore
    return distance(reshaped_matrix1, matrix2)


def calc_vertex_angles(
    vertex: NDArray[np.float32],
    point1: NDArray[np.float32],
    point2: NDArray[np.float32],
    degrees: bool = False,
) -> NDArray[np.float32]:
    """Calculate the angle between three sets of points.

    This function calculates the angle at each `vertex` created by points `point1` and `point2`.

    Args:
        vertex (NDArray[np.float32]): Shape (n, m, 3). The coordinates of the vertex point(s), \
            where the angle is being measured.
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

    normalized_vector1: NDArray[np.float32] = normalize(vector1)
    normalized_vector2: NDArray[np.float32] = normalize(vector2)

    # Dot product of normalized vectors
    dotproduct = dot_product_3d(normalized_vector1, normalized_vector2)

    angle_rad: NDArray[np.float32] = np.arccos(dotproduct)

    if degrees:
        return np.degrees(angle_rad)

    return angle_rad


def calc_vec_angle(vector1: NDArray[np.float32], vector2: NDArray[np.float32]) -> NDArray[np.float32]:
    """Calculate the angle between two vectors in degrees.

    This function computes the angle between pairs of vectors. Each vector in the pair is normalized first,
    and then the dot product is calculated. The angle is obtained using the arccosine of the dot product,
    adjusted by the sign of the dot product. Finally, the angle is converted to degrees.

    Args:
        vector1 (NDArray[np.float32]): The first vector or a set of vectors in each pair.
        vector2 (NDArray[np.float32]): The second vector or a set of vectors in each pair.

    Returns:
        NDArray[np.float32]: The angles between pairs of vectors in degrees.

    ??? example "Example"
        ```py
        v1 = np.array([[1, 0, 0], [0, 1, 0]])
        v2 = np.array([[0, 1, 0], [0, 0, 1]])
        calc_vec_angle(v1, v2)
        array([90., 90.], dtype=float32)
        ```
    """
    vector1 = normalize(vector1)
    vector2 = normalize(vector2)
    dotproduct = dot_product(vector1, vector2)
    raw_angle = np.arccos(np.clip(dotproduct, -1.0, 1.0))
    adjusted_angle = np.sign(dotproduct) * raw_angle
    angle_in_degrees: NDArray[np.float32] = np.degrees(adjusted_angle)
    angle_in_degrees[angle_in_degrees < 0] = 180 + angle_in_degrees[angle_in_degrees < 0]

    return angle_in_degrees


def calc_vec_line_angles(vector: NDArray[np.float32], line_direction: NDArray[np.float32]) -> NDArray[np.float32]:
    """Calculate the angle between a vector and a line direction in degrees.

    This function computes the angle between each pair of vector and line direction.
    The dot product of the vector and the normalized line direction is calculated.
    The raw angle is obtained using the arccosine of the dot product, and it is adjusted by the sign of the dot product.
    If the adjusted angle is less than zero, it is incremented by pi.
    Finally, the angle is converted to degrees.

    Args:
        vector (NDArray[np.float32]): The vector or set of vectors.
        line_direction (NDArray[np.float32]): The direction or set of directions of the lines.

    Returns:
        NDArray[np.float32]: The angles between the vectors and line directions in degrees.

    ??? example "Example"
        ```py
        vec = np.array([[1, 0, 0], [0, 1, 0]])
        line_dir = np.array([[0, 1, 0], [0, 0, 1]])
        calc_vec_line_angles(vec, line_dir)
        array([90., 90.], dtype=float32)
        ```
    """
    normalized_line_direction = normalize(line_direction)
    dotproduct = dot_product(vector, normalized_line_direction)
    raw_angle = np.arccos(dotproduct)
    adjusted_angle = np.sign(dotproduct) * raw_angle

    adjusted_angle = np.where(adjusted_angle < 0, adjusted_angle + np.pi, adjusted_angle)

    result: NDArray[np.float32] = np.degrees(adjusted_angle)
    return result
