"""Choose an available R worker port below Linux's usual ephemeral-port range."""
import random
import socket

# Independent OS randomness keeps simultaneous jobs from choosing adjacent ports.
# It does not read or modify the R random-number state used for simulations.
start = random.SystemRandom().randrange(8000)
for offset in range(8000):
    port = 12000 + (start + offset) % 8000
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
        try:
            probe.bind(("", port))
        except OSError:
            continue
        print(port)
        break
else:
    raise SystemExit("No available R worker port in 12000-19999")
