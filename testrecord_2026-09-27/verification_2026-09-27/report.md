# Integration test results

FINE_ACCEPTED = supplied fine planner accepted the generated path. No global optimality or independent CAD certification is claimed.

| Case | Status | Accepted | Expected match | Seconds |
|---|---|---:|---|---:|
| open_straight | FINE_ACCEPTED | 1 | True | 0.24 |
| open_waypoint | FINE_ACCEPTED | 1 | True | 0.06 |
| blocked_wall | COARSE_NO_PATH | 0 | True | 0.16 |
| blocked_start | FINE_NO_PATH | 0 | True | 0.06 |
| single_obstacle | FINE_ACCEPTED | 1 | None | 0.09 |

Attempts: 5; accepted pairs: 3; test/errors: 0.
