"""
Optimizing the ToF based on some constraints
"""
from typing import Callable, Literal, Tuple


def two_pointer_discrete_optimize(
    target_func: Callable[[int, int], float],
    left_pointer_init: int,
    right_pointer_init: int,
    optima_type: Literal["max", "min"] = "max",
) -> Tuple[int, int, float]:
    # TODO: This function gets stuck in my usecase. There is often times point on the left which produces a local optima
    # Its somehow just a single point. Moving 2 points would solve this issue.
    # TODO: This function is very unoptimized. Optimize it later with larger jumps maybe?
    """
    Use two pointers to optimize the target function by moving one inwards at a time.

    The target function is assumed to always return a float(i.e, not discontinuous) and has a single maximum within
    the pointer limits

    Args:
        target_func (Callable[int, int]): The target function that takes in two integers and returns a float.
        left_pointer_init (int): Starting value(left)
        right_pointer_init (int): Satring value(right)

    Returns:
        List[int, int, float]: The left pointer, right pointer and the optimum value
    """

    def _is_next_better(old, new):
        """
        Extra Note: We want the range to be as large as possible. When given the choise, if the next value is equal, we
        do not consider it as a better choice
        """
        if optima_type == "max":
            return new > old
        elif optima_type == "min":
            return new < old

    left_pointer = left_pointer_init
    right_pointer = right_pointer_init
    optimum_value = target_func(left_pointer, right_pointer)
    while left_pointer < right_pointer:
        next_left = target_func(left_pointer + 1, right_pointer)
        next_right = target_func(left_pointer, right_pointer - 1)
        if _is_next_better(optimum_value, next_left):
            left_pointer += 1
            optimum_value = next_left
        elif _is_next_better(optimum_value, next_right):
            right_pointer -= 1
            optimum_value = next_right
        else:
            break

    return left_pointer, right_pointer, optimum_value
