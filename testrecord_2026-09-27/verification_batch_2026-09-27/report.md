# Integration test results

FINE_ACCEPTED = supplied fine planner accepted the generated path. No global optimality or independent CAD certification is claimed.

| Case | Status | Accepted | Expected match | Seconds |
|---|---|---:|---|---:|
| open_straight | FINE_ACCEPTED | 1 | True | 0.06 |
| open_waypoint | FINE_ACCEPTED | 1 | True | 0.08 |
| blocked_wall | COARSE_NO_PATH | 0 | True | 0.16 |
| blocked_start | FINE_NO_PATH | 0 | True | 0.06 |
| single_obstacle | FINE_ACCEPTED | 1 | None | 0.08 |
| random_0000 | FINE_ACCEPTED | 1 | None | 0.08 |
| random_0001 | FINE_ACCEPTED | 1 | None | 0.08 |
| random_0002 | FINE_ACCEPTED | 1 | None | 0.08 |
| random_0003 | FINE_ACCEPTED | 1 | None | 0.06 |
| random_0004 | FINE_ACCEPTED | 1 | None | 0.06 |
| random_0005 | FINE_ACCEPTED | 1 | None | 0.09 |
| random_0006 | FINE_ACCEPTED | 1 | None | 0.11 |
| random_0007 | FINE_ACCEPTED | 1 | None | 0.08 |
| random_0008 | FINE_ACCEPTED | 1 | None | 0.08 |
| random_0009 | FINE_ACCEPTED | 1 | None | 0.08 |

Attempts: 15; accepted pairs: 13; test/errors: 0.
