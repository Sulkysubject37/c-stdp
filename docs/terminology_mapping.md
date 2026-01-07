# Terminology Mapping for Infectious Disease Symposium

To ensure clinical relevance and avoid neuroscience jargon, the following terminology mapping is strictly enforced in all symposium materials (posters, slides, summary).

| Original Computational Term | Symposium / Clinical Term | Definition / Context |
| :--- | :--- | :--- |
| **Spike** | **Sudden Onset Event** | The specific timepoint where a gene's expression rate of change exceeds a threshold, indicating rapid upregulation. |
| **Spike Train** | **Activation Timeline** | The discrete sequence of onset events for a specific gene over the course of the infection. |
| **STDP (Spike-Timing Dependent Plasticity)** | **Temporal Event Mining** | The core algorithm that infers causality based on the consistent time-delay between activation events. |
| **STDP Edge / Weight** | **Temporal Influence Link** | A directed connection indicating that Gene A consistently activates before Gene B. |
| **Potentiation** | **Consistent Precedence** | The statistical reinforcement of a link when Gene A precedes Gene B. |
| **Depression** | **Sequence Violation** | The weakening of a link when the expected temporal order is reversed. |
| **Neuron / Node** | **Host Factor** | A gene or protein involved in the immune response. |

## Usage Guidelines
- **DO NOT** use "spiking neural network" in the abstract. Use "Event-based temporal inference".
- **DO NOT** discuss "synaptic weights". Discuss "regulatory strength".
