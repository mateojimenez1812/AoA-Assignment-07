"""
Analysis of Algorithms - Assignment 07
Exercise 3: Transform a random array into a heap (bottom-up heap construction)
"""
import random


def sift_down(arr, i, n):
    """Move arr[i] down until the max-heap property holds in the subtree rooted at i."""
    while True:
        largest = i
        left, right = 2 * i + 1, 2 * i + 2
        if left < n and arr[left] > arr[largest]:
            largest = left
        if right < n and arr[right] > arr[largest]:
            largest = right
        if largest == i:
            return
        arr[i], arr[largest] = arr[largest], arr[i]
        i = largest


def build_heap(arr):
    """
    Bottom-up heap construction (in place).
    Starts from the last internal node (n//2 - 1) and sifts each node down.
    Running time: Theta(n).
    """
    n = len(arr)
    for i in range(n // 2 - 1, -1, -1):
        sift_down(arr, i, n)
    return arr


def is_max_heap(arr):
    n = len(arr)
    for i in range(n):
        for child in (2 * i + 1, 2 * i + 2):
            if child < n and arr[i] < arr[child]:
                return False
    return True


def print_tree(arr):
    """Print the heap level by level."""
    level, start = 0, 0
    while start < len(arr):
        size = 2 ** level
        print(f"   level {level}: {arr[start:start + size]}")
        start += size
        level += 1


if __name__ == "__main__":
    print("EXERCISE 3: Transform a random array into a heap")
    random.seed(3)
    arr = [random.randint(1, 99) for _ in range(12)]
    print("Random array :", arr)
    print("Is max-heap? :", is_max_heap(arr))
    build_heap(arr)
    print("\nAfter build_heap:", arr)
    print("Is max-heap? :", is_max_heap(arr))
    print("\nHeap by levels:")
    print_tree(arr)

    ok = 0
    for _ in range(1000):
        a = [random.randint(-1000, 1000) for _ in range(random.randint(0, 300))]
        original = sorted(a)
        build_heap(a)
        if is_max_heap(a) and sorted(a) == original:
            ok += 1
    print(f"\n{ok}/1000 random arrays transformed into valid heaps (same elements kept)")
