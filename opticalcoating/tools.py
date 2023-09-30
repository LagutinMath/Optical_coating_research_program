import functools, time

def timer(iters=1):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            total = 0
            for i in range(iters):
                start = time.perf_counter()
                value = func(*args, **kwargs)
                end = time.perf_counter()
                total += end - start
            print(f'Среднее время выполнения {func.__name__}: {total/iters:.2} сек.')
            return value
        return wrapper
    return decorator