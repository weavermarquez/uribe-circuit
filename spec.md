# Uribe Oscillator Network Simulation — Comprehensive Specification

## 1. Motivation and Purpose

This project simulates the oscillating network described in Ricardo Uribe's *Tractatus Paradoxico-Philosophicus*. Uribe, a collaborator with Maturana and Varela on autopoiesis, physically built a network of coupled oscillators on a 10×10 grid to demonstrate how synchronized, coherent behavior **emerges** from simple local interactions — and how this emergence serves as a metaphor for living systems, societies, language, and the relationship between logical and paradoxical reasoning.

**This is not a web app or a game.** It is a simulation of a physical system for scientific demonstration. The primary goal is behavioral fidelity to Uribe's described (and physically constructed) circuit. Performance is secondary — modern hardware is more than sufficient. Visual appeal is tertiary — correctness comes first.

### What Uribe Built

A 10×10 grid of oscillator nodes, each consisting of a logical inverter in a self-referential feedback loop (a NOT gate whose output feeds back to its input through a propagation delay). Each node independently oscillates between ON (red/light) and OFF (white/dark). When a node transitions to OFF, it sends an **excitatory signal** to its four Von Neumann neighbors (up, down, left, right), encouraging them to turn ON. The grid wraps toroidally (edges connect to opposite edges), and the wrap-around always connects even rows/columns to odd ones, maintaining bipartite structure.

From random initial conditions, the network **self-synchronizes** into a checkerboard pattern: alternating ON/OFF states between adjacent nodes, with the two sublattices anti-phase locked. This emergent synchronization is the core phenomenon to replicate.

### Why It Matters

Uribe uses this network to demonstrate:

1. **Dynamic stability vs. static stability**: The checkerboard is not a fixed state but an ongoing *oscillation* maintained by continuous activity. Stop the activity and the pattern vanishes. This is his metaphor for living systems (autopoiesis) vs. non-living things.
2. **Observer-dependence without observer-influence**: An observer can trace paths, follow diagonals, select groups — but these choices don't affect the network. Yet without an observer, the patterns "make no sense." The network just oscillates.
3. **Robustness and fragility**: Local damage (a few dead nodes) is tolerated. Major damage causes the coherent oscillation to disintegrate. This maps to his claims about societies, language, and hierarchies.
4. **Paradoxical and logical reasoning**: The oscillation itself embodies the paradox (the NOT gate that is both true AND false when considered paradoxically, but alternates true OR false when considered logically in time).

---

## 2. Reference Material

### Primary Source

Ricardo Uribe, *Tractatus Paradoxico-Philosophicus*, particularly:
- The "Metaphor" section (network of oscillators)
- The "Introduction" section (from logics to paradoxes, paradoxical toroids, the DPDT switch example)

Available at: `http://bcl.ece.illinois.edu/Uribe/`

### Key Quotes Governing Behavior

> "A basic oscillator consists of a logical inverter within a paradoxical loop (a loop with a 'twist' such as the Möbius band). As with a paradox the oscillation develops in a time and a space resulting from the observer contemplating the oscillator from a logical perspective."

> "When a node turns off its light, it sends an excitatory signal to its four neighbors (up, down, right and left, including a wrap-around at the borders of the matrix) so that their lights turn on."

> "This activity coming from the nodes synchronizes the whole network so that it oscillates turning on and off alternate lights (vertically and horizontally)"

> "The choices of the observer do not affect the oscillation of the network but without the network's activity the observer cannot choose."

> "may sustain some local damage (some nodes stop oscillating, some excitatory or inhibitory signals blocked, etc.) but the network as a whole continues its dance... With major damage the activity of the network may collapse."

> "The wrap-around, structural or not, operates from an even row (column) to an odd one or from an odd to an even one."

> "it can interact with another network or with itself (e.g., using a mirror)"

> "Inhibiting the activity of some lights, an external stimulus (e.g., light) may lead to the partial or complete isolation of other lights. These will continue to oscillate by themselves, but the activity of the network will disintegrate and cease to oscillate as an entity."

### Uribe's Physical Construction

Uribe states: "(I built this network and it operates as described)." He also built a multi-dimensional version and the DPDT switch circuit from the Introduction. The exact schematic is not provided, but the description implies:

- Logic-level oscillators (NOT gate + feedback delay), not analog/sinusoidal
- Human-visible oscillation rates (Hz range), suggesting either slow analog components (RC, 555 timers) or frequency-divided digital logic
- Physical wire connections between adjacent nodes on a PCB or breadboard
- Toroidal wrap-around via wiring

---

## 3. Core Behavioral Invariants

These are the **acceptance criteria**. If the simulation does not produce these behaviors, it is incorrect. Each invariant maps to a specific passage in Uribe and must have corresponding automated tests.

### INV-1: Autonomous Oscillation
Each node oscillates independently. Remove all coupling (disconnect all edges) and every node still flips ON→OFF→ON→OFF forever at its intrinsic frequency. The oscillation is NOT driven by neighbors — it is intrinsic (the "paradoxical loop" / self-referential NOT gate).

**Test**: Isolated node with no neighbors oscillates at period 2τ indefinitely. No drift, no decay.

### INV-2: Excitatory Coupling on OFF-Transition Only
When a node transitions from ON to OFF, it sends an excitatory signal to all its neighbors. This signal encourages neighbors to turn ON. There is NO signal sent on the ON transition. There is NO inhibitory coupling between nodes — the only "inhibition" is each node's internal self-flip.

**Test**: Node A goes OFF → neighbors receive pulse. Node A goes ON → no pulse emitted. Verify by inspecting event log.

### INV-3: Anti-Phase Checkerboard Emergence
Starting from random initial conditions (random states, random phase offsets), the network converges to a checkerboard pattern where adjacent nodes are in opposite states. The two sublattices (even-parity and odd-parity nodes, where parity = (x+y) mod 2) become anti-phase locked.

**Test**: Start random, compute neighbor-disagreement ratio (fraction of adjacent pairs in opposite states) at each half-period. Verify convergence toward 1.0 within a reasonable number of cycles (~5-20 depending on parameters).

### INV-4: Observer Independence
Any visualization, path-tracing, group-selection, highlighting, or analysis is PURELY read-only with respect to the simulation. No observer action can modify the event log or network state.

**Test**: Same seed and parameters produce bit-identical event logs regardless of renderer state, observer tool activity, or whether a renderer is attached at all.

### INV-5: Self-Sustaining Activity via Internal Feedback
There is no global clock and no external driver. The network's continued oscillation is sustained entirely by (a) each node's intrinsic oscillation and (b) local coupling between neighbors. Activity anywhere leads to further activity elsewhere through feedback loops.

**Test**: Network runs indefinitely without external input. Removing coupling → nodes desync but each continues to oscillate. Freezing all nodes (disabling intrinsic oscillation) → all activity stops.

### INV-6: Local Damage Tolerance, Global Damage Fragility
Killing a small number of nodes (setting them to permanently OFF, no signals) does not destroy the checkerboard pattern — the network adapts around the gaps. Killing a large fraction of nodes causes the coherent oscillation to disintegrate.

**Test**: Kill 5% of nodes → coherence metric recovers above 0.9. Kill 50% → coherence falls below 0.6. Determine approximate percolation threshold.

### INV-7: Toroidal Wrap-Around Maintains Bipartiteness
The wrap-around boundary condition connects even rows/columns to odd rows/columns, preserving the bipartite graph structure necessary for checkerboard formation.

**Test**: For a 10×10 grid, verify that node (9, y) is neighbors with node (0, y), and that these always have opposite sublattice parity. Verify the full graph is bipartite.

### INV-8: Observable Diagonal Wavefronts
In steady state, ON cells form diagonal lines that appear to translate by one position per half-period. These are a stroboscopic consequence of the checkerboard alternation, not a separate propagation mechanism.

**Test**: In steady state, extract connected diagonal runs of ON cells at consecutive half-periods. Verify spatial shift of +1 or -1 per step.

### INV-9: External Perturbation Modifies Dynamics
External stimuli (forcing nodes ON/OFF, injecting signals) alter the network's attractor landscape. Some neighbors may not turn on "when expected," changing the observable patterns.

**Test**: Inject perturbation → observe local desynchronization in coherence metric → observe recovery (or persistent alteration).

### INV-10: Network Self-Interaction
Output signals from boundary nodes can be routed back as input to the same or other networks. "Mirror mode": right-edge output → left-edge input.

**Test**: Enable mirror routing → observe modified steady-state dynamics compared to standard toroidal wrap.

### INV-11: Isolation via Inhibition Fragments the Entity
Blocking all coupling across a boundary (e.g., a vertical line of edges) isolates regions. Each region continues to oscillate internally (INV-1) but decouples from the whole. The network ceases to function "as an entity."

**Test**: Block all edges crossing x=5 → left half and right half each independently converge to their own checkerboard (possibly different phase), but desynchronize from each other.

---

## 4. System Model

### 4.1 Node Model — Binary Logical Oscillator

Each node is a NOT gate with self-referential feedback. The oscillation arises from propagation delay, not from a continuous phase variable.

```
Node {
  id:           NodeId          // (x, y) coordinates or generic ID
  state:        ON | OFF        // binary, no intermediate values
  tau:          float           // propagation delay (time for NOT gate to flip)
  t_last_flip:  float           // sim-time of most recent state change
  generation:   int             // incremented on every reschedule (for event cancellation)
  alive:        boolean         // false = permanently dead, no oscillation or signals
}
```

**In isolation**, a node oscillates as a square wave:
- ON for duration τ, then flips to OFF
- OFF for duration τ, then flips to ON
- Period = 2τ

The intrinsic oscillation period comes from `state(t + τ) = NOT state(t)`. This is the "paradoxical loop."

**Heterogeneity**: Each node's τ is drawn from `N(τ₀, σ_τ²)` at initialization, representing manufacturing variation in the physical circuit. This heterogeneity is essential — without it, a perfectly regular grid with identical delays would sync trivially or not at all, with no observable *process* of emergence.

### 4.2 Signal Model — Excitatory Pulse

When node `i` transitions **ON → OFF**, it emits one excitatory pulse to each neighbor. The pulse arrives after wire delay `δ`.

**On pulse arrival at node `j`:**

| Node j state | Effect |
|---|---|
| OFF | **Strong coupling**: Immediately flip j to ON. Reset j's internal timer to τ_j (full ON duration from now). Cancel j's pending SELF_FLIP, schedule new one. |
| OFF | **Weak coupling**: Reduce j's remaining OFF-time by factor κ. `t_remaining *= (1 - κ)`. Reschedule SELF_FLIP. |
| ON | No effect. (Node is already in the desired state. Optionally: reset/extend ON timer.) |

Both coupling modes should be implemented. Strong coupling (κ = 1.0) is the digital interpretation. Weak coupling (κ < 1.0) is the analog approximation. Uribe's physical circuit behavior is uncertain — parameterize and let the user compare.

**Why there are no race conditions**: The signal space is asymmetric. Excitatory pulses can only push nodes toward ON. The only force pushing toward OFF is the node's internal self-flip. These cannot conflict: a pulse arriving at an ON node is a no-op; a pulse arriving at an OFF node accelerates its return to ON. With continuous time, heterogeneous τ, and nonzero δ, exact simultaneous events are measure-zero.

**Refractory period**: There is no separate refractory parameter. The propagation delay τ IS the refractory period. After a node flips to ON (whether from self-flip or from pulse), it cannot flip back to OFF for τ — that's how long the NOT gate takes to propagate. A pulse arriving at an ON node simply has no effect.

### 4.3 Coupling Variants to Explore

Uribe's text does not specify the exact coupling mechanism at the circuit level. Three physical interpretations exist:

1. **DC-coupled (level)**: Neighbor's input is held continuously by sender's output. Risk: latching (neighbors hold each other in complementary states indefinitely, suppressing oscillation). This is probably NOT what Uribe built, since the network keeps oscillating.

2. **AC-coupled (pulse/edge)**: Only transitions propagate. This matches "sends an excitatory signal" language. The default model.

3. **Edge-triggered**: Only the transition edge matters, with fixed pulse width. Similar to AC-coupled but more digital.

All three should be implementable by varying the coupling handler. Start with AC-coupled (option 2). Parameterize.

### 4.4 Topology

```
Graph {
  nodes:      Map<NodeId, Node>
  adjacency:  Map<NodeId, Set<NodeId>>   // undirected, mutable at runtime
}
```

**Topology builders:**
- `toroidalGrid(rows, cols)`: 2D torus with Von Neumann neighborhood. Wrap-around connects even↔odd.
- `hypercubicTorus(dimensions, size)`: n-D torus. Node has 2n neighbors. Generalizes the 2D case.
- `pair()`: 2-node test graph.
- `ring(n)`: cycle graph (1D torus).
- `custom(adjacencyList)`: arbitrary graph for irregular topologies.

**Runtime mutability**: The graph must support `removeEdge(a, b)`, `killNode(id)`, and `addEdge(a, b)` for perturbation experiments (INV-6, INV-9, INV-10, INV-11). These operations take effect immediately in sim-time.

**Bipartiteness**: A toroidal grid is bipartite IFF all dimensions are even. Odd dimensions create frustrated cycles where a clean checkerboard is impossible. This is not a bug — it is an interesting experimental condition. The UI should indicate whether the current topology is bipartite, but not enforce even dimensions.

### 4.5 Time Model — Simulation Time Only

The simulation operates exclusively in **sim-time**. There is no wall-clock coupling in the engine. The engine produces an ordered event log; a renderer (if attached) maps sim-time to wall-clock time for playback.

**Event types:**

```
type SimEventType = 'SELF_FLIP' | 'PULSE_EMIT' | 'PULSE_ARRIVE' | 'NODE_KILL' | 'EDGE_BLOCK' | 'EXTERNAL_FORCE'

type SimEvent = {
  t:          number              // sim-time
  id:         number              // global monotonic event counter (tiebreaker)
  type:       SimEventType
  target:     NodeId
  source?:    NodeId              // for PULSE_ARRIVE: which node sent it
  generation: number              // must match target node's generation to be valid
  metadata?:  Record<string, any> // coupling mode, forced state, etc.
}
```

**Event priority**: Primary sort by `t` (ascending). Secondary sort by `id` (ascending, FIFO within same timestamp). This provides deterministic ordering.

**Event cancellation via lazy deletion**: Rather than searching the heap to remove cancelled events:
- Each node holds a `generation` counter.
- Each event stores the `generation` at creation time.
- When a node is rescheduled (e.g., pulse forces it ON, resetting its timer), increment its generation.
- When dequeuing an event, if `event.generation ≠ node.generation`, discard it silently.

This is O(1) cancellation with O(log N) overhead on dequeue (discarding stale events).

### 4.6 Engine Loop

```
function step():
  event = queue.dequeue()
  if event is stale (generation mismatch): discard, return null
  advance sim_clock to event.t

  switch event.type:
    case SELF_FLIP:
      flip node state
      node.t_last_flip = sim_clock
      if new state is OFF:
        for each neighbor in adjacency[node.id]:
          enqueue PULSE_ARRIVE at neighbor, time = sim_clock + δ
      increment node.generation
      enqueue next SELF_FLIP at time = sim_clock + node.tau, with new generation

    case PULSE_ARRIVE:
      if target node is dead: discard
      if target node is ON: discard (no-op)
      if target node is OFF:
        if strong coupling:
          flip target to ON
          target.t_last_flip = sim_clock
          increment target.generation
          enqueue SELF_FLIP at time = sim_clock + target.tau, with new generation
        if weak coupling:
          reduce remaining OFF time by factor κ
          reschedule pending SELF_FLIP accordingly
          increment target.generation

    case NODE_KILL:
      set node.alive = false
      increment generation (invalidates all pending events)

    case EDGE_BLOCK:
      remove edge from adjacency

    case EXTERNAL_FORCE:
      force node to specified state
      reset timer accordingly
      increment generation

  append event to event_log
  return event

function run(t_end):
  while queue.peek().t <= t_end:
    step()
  return event_log
```

### 4.7 Parameters

| Parameter | Symbol | Default | Description |
|---|---|---|---|
| Base propagation delay | τ₀ | 1.0 | Mean time for NOT gate to flip (arbitrary time units) |
| Delay heterogeneity | σ_τ | 0.01 (1%) | Std dev of τ across nodes |
| Wire delay | δ | 0.02 | Time for pulse to travel between adjacent nodes |
| Coupling strength | κ | 1.0 | 1.0 = strong (instant force-ON); <1.0 = weak (proportional) |
| Grid rows | rows | 10 | |
| Grid columns | cols | 10 | |
| Dimensions | n_dim | 2 | For hypercubic torus |
| Random seed | seed | any integer | For reproducible runs |

**Time units are arbitrary.** τ₀ = 1.0 means one time unit per half-period, full period = 2.0. The renderer maps sim-time to wall-clock time via a speed multiplier.

### 4.8 Initialization

1. Build graph topology from parameters.
2. Seed PRNG.
3. For each node:
   - `tau = τ₀ + PRNG.gaussian() * σ_τ` (clamp to positive)
   - `state = PRNG.random() < 0.5 ? ON : OFF` (random initial state)
   - `t_last_flip = 0`
   - `generation = 0`
   - `alive = true`
   - Schedule initial SELF_FLIP at time `PRNG.random() * tau` (random phase offset within first half-period)
4. Sim clock starts at 0.

### 4.9 Observables and Metrics

**Coherence (anti-phase order parameter)**:
```
coherence = (number of adjacent pairs in opposite states) / (total number of adjacent pairs)
```
- 1.0 = perfect checkerboard
- 0.5 = random
- 0.0 = perfect in-phase (all same state) — should not occur with excitatory coupling

Sample at regular sim-time intervals (e.g., every τ₀) by inspecting current node states.

**Time-to-sync**: Number of full periods (2τ₀) until coherence first exceeds 0.95.

**Node state history**: For any selected node, the sequence of (t, state) pairs from the event log.

**Diagonal detection**: In steady state, identify connected diagonal runs of ON cells and track their spatial shift between consecutive half-periods.

---

## 5. Architecture

### 5.1 Separation of Concerns

```
SimEngine (pure, no side effects, no DOM, no rendering)
  ├── EventQueue         (min-heap with lazy deletion via generation counters)
  ├── Graph              (generic adjacency list, mutable)
  ├── Node[]             (state, tau, generation, t_last_flip)
  ├── PRNG               (seeded, deterministic)
  ├── step() → SimEvent | null
  ├── run(t_end) → SimEvent[]
  └── getState() → NetworkState (snapshot)

EventLog (append-only, indexed)
  ├── events: SimEvent[]
  ├── stateAt(t) → NetworkState      // binary search + replay from last snapshot
  ├── coherenceAt(t) → number
  └── toJSON() / fromJSON()

Renderer (reads EventLog or live SimEngine state, NEVER writes to sim)
  ├── GridCanvas         (node colors, phase brightness, transition indicators)
  ├── PlaybackController (play, pause, speed, scrub, step-forward, step-backward)
  ├── MetricsDisplay     (coherence number + line chart)
  └── future: ObserverTools, PerturbationUI, 3D views

Controller (bridges sim ↔ renderer ↔ UI)
  ├── TimeGovernor       (maps sim-time to wall-clock via speed multiplier)
  ├── ParameterPanel     (sliders, toggles, seed input)
  └── PerturbationDispatcher (converts UI actions to sim events)
```

### 5.2 Key Constraints

- **SimEngine must be runnable headlessly.** No browser APIs. No DOM. No canvas. Pure TypeScript logic. Testable with `bun test`.
- **EventLog is the source of truth.** Renderer reads from it. SimEngine writes to it. Nothing else touches it.
- **Observer tools are read-only overlays.** They query the EventLog but never modify it or the SimEngine (INV-4).
- **Perturbations go through the engine.** UI actions (kill node, block edge) are converted to sim events and enqueued. They don't directly mutate state.

### 5.3 Time Governor (Renderer Only)

The renderer maps sim-time to wall-clock time:
```
each render frame:
  wall_dt = performance.now() - last_frame_time
  sim_advance = wall_dt * speed_multiplier / 1000  // convert ms to sim-time-units
  engine.run(sim_clock + sim_advance)
  render engine.getState()
```

**Speed multiplier**: 1.0 = real-time (if τ₀=1.0, one sim-second = one wall-second). 10.0 = 10x fast-forward. 0.1 = slow-mo. Adjustable via UI slider.

**Sub-frame transitions**: Between two render frames (16.67ms at 60fps), multiple state transitions can occur on a single node. The renderer shows **state at the current sim-time**, not intermediate states. Optional: transition indicator (brief glow) for nodes that changed since last frame.

### 5.4 Determinism and Reproducibility

- All stochastic elements use the seeded PRNG.
- The seed is displayed in the UI and can be set manually.
- Same seed + same parameters = identical event log, guaranteed.
- The event log can be serialized to JSON and replayed.

---

## 6. Technology Stack

### 6.1 Runtime and Language

- **Bun** as runtime (fast, TypeScript-native, built-in test runner)
- **TypeScript** throughout (strict mode)

### 6.2 Dependencies (Minimal)

| Dependency | Purpose | Why not hand-roll |
|---|---|---|
| `flatqueue` or `tinyqueue` | Min-heap priority queue for event scheduling | Battle-tested, typed-array-backed, saves 30 min vs hand-rolling a correct heap |
| `seedrandom` | Seeded PRNG for reproducibility | Mature, drop-in `Math.random()` replacement, multiple algorithms |

### 6.3 What We Are NOT Using

| Category | Rejected Options | Rationale |
|---|---|---|
| Circuit simulators | SPICE, ngspice, Xyce | Solve continuous differential equations (Kirchhoff's laws). Our system is fundamentally discrete-event. SPICE's timestep-based solver is the wrong paradigm. Integration pain exceeds benefit. |
| Digital HDL simulators | Verilator, Icarus Verilog, digital.js | Model logic gates with propagation delays (close to our needs!). But HDL is a detour from the scientific questions; coupling between "circuits" isn't natural in HDL; these tools optimize for verification not visualization. |
| Mixed-signal simulators | — | Worst of both worlds for our purposes. |
| Graph libraries | ngraph.graph, graphology | Our topology needs (n-D torus + mutable edges) are simple enough that a hand-rolled adjacency list is cleaner and avoids dependency overhead. |
| Heavy frameworks | Rama, Kafka, etc. | 100 nodes generating ~10k events/sim-second. A min-heap in TypeScript handles millions of ops/second. Infrastructure overhead vastly exceeds compute needs. |
| D3 / chart libraries | d3, chart.js | Overkill for initial metrics display. A simple canvas line chart suffices. Add later if needed. |

### 6.4 Rendering

- **Canvas 2D** for the grid visualization (sufficient for 2D grids up to ~100×100)
- Plain HTML/CSS for UI controls
- Upgrade path to WebGL if/when 3D visualization or larger grids are needed (PR 6+)
- No React. No framework. Vanilla DOM manipulation for controls. The simulation is the complex part, not the UI.

### 6.5 Testing

- **bun:test** — zero-config, built into runtime
- SimEngine tests run headlessly (no browser, no DOM)
- Each invariant (INV-1 through INV-11) has at least one corresponding test
- Tests are deterministic (seeded PRNG)

### 6.6 Build

- **No build step initially.** Bun runs TypeScript directly.
- Add Vite or similar only when browser builds are needed (PR 2+).

---

## 7. PR Plan — Vertical Slices

Each PR delivers a testable, demonstrable increment. The minimal viable renderer is included from PR 1 so that every slice is visually verifiable.

### PR 1: Core Engine + Validation + Terminal Renderer

**Scope**: The simulation engine, topology builders, automated tests for core invariants, and a minimal terminal-based renderer for visual verification.

**SimEngine:**
- `Node` type: id, state, tau, t_last_flip, generation, alive
- `Graph`: mutable adjacency list (Map<NodeId, Set<NodeId>>), addNode, removeNode, addEdge, removeEdge
- `EventQueue`: min-heap (via flatqueue), lazy deletion via generation counters, global monotonic event ID for tiebreaking
- `SimEngine.step()` → processes one event, appends to log, returns event or null (if stale)
- `SimEngine.run(t_end)` → processes all events up to t_end, returns event log slice
- `SimEngine.getState()` → current NetworkState snapshot (all node states)
- `SimEngine.getCoherence()` → current anti-phase coherence metric
- Seeded PRNG (seedrandom) for all stochastic initialization
- Strong coupling mode (κ = 1.0): pulse forces OFF→ON immediately

**Topology builders:**
- `makePair()` → 2-node graph (minimal test case)
- `makeRing(n)` → cycle graph (1D torus)
- `makeToroidalGrid(rows, cols)` → 2D torus with Von Neumann neighborhood

**Tests (bun:test):**
- **INV-1**: Isolated node (no neighbors) oscillates at period 2τ for 100+ cycles. Verify exact flip times.
- **INV-2**: 2-node pair: verify pulse emitted only on OFF transition. Inspect event log for absence of pulse on ON transition.
- **INV-3**: 10×10 grid from random seed: coherence exceeds 0.95 within 20 full periods.
- **INV-4**: Same seed → identical event log across two runs.
- **INV-5**: Remove all edges mid-run → nodes continue oscillating individually. Freeze all nodes → no further events.
- **INV-7**: Verify bipartiteness of 10×10 toroidal grid (all edges cross sublattice parity).
- 2-node anti-phase lock: verify nodes reach opposite states within 5 cycles.
- 4-node ring: verify alternating pattern on smallest bipartite cycle.

**Terminal renderer:**
- Script that replays event log as ASCII frames to stdout.
- Grid displayed as characters: `█` for ON, `·` for OFF.
- Print frame at regular sim-time intervals (e.g., every τ₀/2).
- Show sim-time and coherence metric per frame.
- Optional: write event log to JSON file for external analysis.

**Acceptance criteria:**
- All tests pass.
- Terminal renderer shows visible convergence from random start to checkerboard on 10×10 grid.

### PR 2: Canvas Renderer + Playback Controls

**Scope**: Browser-based visualization that reads from the SimEngine (live) or EventLog (replay).

**Renderer:**
- HTML canvas 2D, sized to grid
- Node cells colored by state: ON = red, OFF = white
- Phase brightness: cell opacity modulated by `(sim_now - t_last_flip) / tau` — newly flipped nodes are bright, about-to-flip nodes are dim. This provides visual phase information without breaking the binary model.
- Transition indicator: brief border flash on nodes that changed state since last render frame, decaying over ~100ms wall-time.

**Playback:**
- TimeGovernor: maps sim-time to wall-clock via requestAnimationFrame + speed multiplier
- Controls: play/pause button, speed slider (0.1x to 100x), step-forward (advance one event), reset
- Scrub bar: drag to any point in sim-time (requires EventLog.stateAt(t) reconstruction)

**Metrics display:**
- Current coherence value as text
- Mini line chart of coherence over sim-time (canvas-drawn, no library)

**Acceptance criteria:**
- Visually matches terminal renderer output for same seed.
- Smooth playback at 60fps for 10×10 grid at 1x speed.
- Scrubbing works without artifacts.

### PR 3: Parameter Controls + Coupling Modes

**Scope**: UI for adjusting simulation parameters and support for weak coupling.

**Parameters UI:**
- Sliders: τ₀, σ_τ, δ, κ, rows, cols
- Seed input: text field, display current seed
- Reset button: re-initialize with current parameters
- All parameter changes require reset (no hot-reloading mid-sim for determinism)

**Weak coupling mode:**
- When κ < 1.0, pulse shortens remaining OFF time proportionally instead of forcing instant ON.
- Implementation: on PULSE_ARRIVE at OFF node, compute `t_remaining = scheduled_flip_time - sim_now`, set new flip time to `sim_now + t_remaining * (1 - κ)`, increment generation, reschedule.

**Acceptance criteria:**
- κ = 1.0 matches PR 1/2 behavior exactly.
- κ = 0.1 shows visibly slower convergence.
- σ_τ = 0.1 (10%) shows messier transients than σ_τ = 0.01 (1%).
- Same seed + same params = same result across runs.

### PR 4: Perturbations + Damage

**Scope**: Interactive tools for testing robustness invariants.

**Perturbation tools:**
- Click node to kill (toggle alive/dead). Dead nodes: gray, no oscillation, no signals.
- Shift-click edge to block. Blocked edges: visually indicated (dashed line or gap).
- Right-click node: force ON or force OFF (external stimulus, one-shot).
- "Isolate region" tool: click-drag to select rectangular region, block all edges crossing the selection boundary (INV-11 test).

**Perturbation implementation:**
- All perturbations are enqueued as sim events (NODE_KILL, EDGE_BLOCK, EXTERNAL_FORCE) at current sim-time.
- They flow through the normal event processing pipeline.
- Event log records all perturbations for replay.

**Tests:**
- **INV-6**: Kill 5 random nodes → coherence recovers above 0.9 within 20 periods. Kill 50 random nodes → coherence stays below 0.6.
- **INV-9**: Force one node to fixed ON → observe local disruption, measure recovery time.
- **INV-11**: Block vertical line of edges → two halves desync (coherence computed separately for each half remains high; cross-boundary coherence drops).

**Acceptance criteria:**
- Perturbations visually disrupt then (for local damage) allow recovery.
- Major damage visibly disintegrates coherent oscillation.
- Perturbations are recorded in event log and survive replay.

### PR 5: Observer Tools

**Scope**: Interactive read-only tools for exploring the network as Uribe describes.

**Node inspector:**
- Click node → highlight it and its neighbors.
- Show time-series of its state (square wave plot) in a side panel.
- Show current τ, frequency, time since last flip.

**Path tracer (Uribe's "follow one light"):**
- Click a node to start. When it goes OFF, its four neighbors highlight as choices.
- Click one of the four → trace continues. Path drawn on grid as colored line.
- Multiple concurrent paths in different colors.
- Paths are purely visual overlays, zero effect on sim (INV-4).

**Group tracker:**
- Select 2+ nodes → highlight as a group.
- Show collective state pattern over time.
- Detect "seesaw" and "windmill" patterns as Uribe describes.

**Diagonal detector:**
- Auto-detect diagonal runs of ON cells in steady state.
- Overlay arrows showing direction of apparent diagonal motion.
- Selectable: track specific diagonals across time.

**Acceptance criteria:**
- All observer tools function without affecting sim output (verified by determinism test).
- Path tracer allows the "choose among four neighbors" interaction Uribe describes.

### PR 6: Advanced Topology + Self-Interaction

**Scope**: Multi-dimensional grids, mirror mode, irregular topologies.

**n-D hypercubic torus:**
- `makeHypercubicTorus(dimensions: number[], size: number)` → Graph
- 3D: 6 neighbors per node. 4D: 8 neighbors. n-D: 2n neighbors.
- Renderer: for 3D, isometric projection or layered 2D slices. For 4D+, projections or slice views.

**Mirror mode (INV-10):**
- Configurable signal routing: output from one boundary can be fed as input to another boundary (same network or different network instance).
- Default mirror: right-edge OFF-transition pulses are delivered to left-edge nodes (and vice versa), INSTEAD of toroidal wrap-around. This is a different topology, not an addition to toroidal wrap.

**Irregular topologies:**
- Random graph: Erdős–Rényi with configurable edge probability.
- Scale-free: Barabási–Albert preferential attachment.
- Custom: load adjacency list from JSON.
- Odd-dimension grids: observe frustrated dynamics (not a bug).

**Acceptance criteria:**
- 3D grid with even dimensions produces 3D checkerboard pattern.
- Mirror mode produces visibly different steady-state from toroidal mode.
- Odd-dimension grid shows persistent frustration (coherence plateaus below 1.0).

---

## 8. Open Questions and Future Directions

### 8.1 Unsettled Design Decisions

- **DC vs. AC vs. edge-triggered coupling**: Uribe's text doesn't specify. All three should eventually be testable. AC-coupled (pulse on transition) is the default.
- **Pulse-on-ON-transition too?**: Current model only emits on OFF→. A variant where ON→ also emits (inhibitory? or second excitatory?) might match some circuit implementations. Parameterize.
- **Multiple simultaneous pulses**: If a node receives pulses from 2+ neighbors within a very short window, should they compound? Current model: first pulse flips to ON, subsequent ones are no-ops (node is already ON). This seems correct for strong coupling. For weak coupling, multiple pulses could additively reduce OFF time.
- **Float precision for long runs**: After 10^6 cycles, accumulated float error in event scheduling could cause drift. Mitigation: compute `t_next = t_base + n * τ` (epoch-relative) rather than `t_next = t_now + τ` (incremental). Not critical for initial PRs.

### 8.2 Future Extensions

- **Sound**: Map node states to audio — each node as a tone, network as a chord. Uribe mentions "dance" repeatedly.
- **Multi-network interaction**: Two networks coupled via boundary signals.
- **Evolutionary dynamics**: Let τ values or coupling strengths mutate over time.
- **Recording and export**: Save event logs, export animation frames, generate publication-quality figures.
- **Comparison mode**: Run two simulations side-by-side with different parameters.

### 8.3 Relation to Autopoiesis

This simulation is a concrete instantiation of Maturana and Varela's autopoietic theory (Uribe was their collaborator). The network is a "closed network of processes" that "produces and maintains itself." The checkerboard pattern is not imposed — it is produced by the network's own activity and maintained by the same activity. Destroying the activity destroys the pattern. This is the dynamic stability at "the heart of all living organisms" that Uribe describes, as opposed to the static stability of non-living things (the "logical toroid" vs. the "paradoxical toroid").

The simulation should make this viscerally obvious: you watch order emerge from chaos, sustained by nothing but local interactions, fragile to sufficient damage, robust to minor perturbation. The observer can trace paths through it but cannot alter it. This is the pedagogical purpose.
