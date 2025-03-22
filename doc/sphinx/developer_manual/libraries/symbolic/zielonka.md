# Zielonka

The standard Zielonka solving algorithm, defined in Algorithm requires that every vertex has an outgoing edge.
If this is the case then the graph is called *total*.
We can achieve this by extending every disjunctive PBES equation with $\textsf{true} \land X_\textsf{false}$ where $X_\textsf{false}$ is defined as $\mu X_\textsf{false} = \textsf{true} \land X_\textsf{false}$ and similarly for the conjunctive PBES equations with $X_\textsf{true}$.
However, adding these unnecessary transitions can be costly.

Therefore, if we perform the deadlock detection we can avoid extending the PBES and obtain a total graph by performing the preprocessing step described in Algorithm~\ref{alg:preprocess}.
Every disjunctive vertex that is a deadlock is won by player odd (previously indicated by a transition to $X_\textsf{false}$) and every conjunctive vertex that is a deadlock by player even.
If we compute the attractors to these won vertices then the resulting graph is total and can be used in the Zielonka algorithm as follows $\textsf{zielonka}(\textsf{preprocess}(V, D, \emptyset, \emptyset))$.


```{imgmath}
\usepackage{algorithm}

\begin{algorithmic}[1]
\Function{Preprocess}{$V, D, W_0, W_1$}
	\State {$W'_0 \gets W_0 \cup (D \cap V_1)$}	
	\State {$W'_1 \gets W_1 \cup (D \cap V_0)$}
	
	\State {$W'_0 \gets \attr{0}{W'_0, V}$}
	\State {$W'_1 \gets \attr{1}{W'_1, V}$}
	
	\State \Return $V \setminus (W'_0 \cup W'_1)$
\EndFunction
\end{algorithmic}
```

```{imgmath}
\usepackage{algorithm}

\begin{algorithmic}[1]
    \Function{Zielonka}{$V$}
        \If {$V = \emptyset$}
            \State \Return {$\emptyset, \emptyset$}
        \EndIf
        \State $m:=\min (\{r(v)\mid v\in V\})$ 
        \State $\alpha := m \bmod 2$
        \State $U:=\{v\in V\mid r(v)=m\}$ 
        \State $A := \attr{\alpha}{U,V}$
        \State $W_{0}^{\prime },W_{1}^{\prime}:= \textsc{Zielonka}(V\setminus A)$
            \If {$W_{1-\alpha }^{\prime }= \emptyset$}
        \State $W_{\alpha },W_{1-\alpha }:=A\cup W_{\alpha }^{\prime},\emptyset$ 
        \Else
            \State $B := \attr{1-\alpha}{W_{1-\alpha }^{\prime }, V}$
            \State $W_{0},W_{1}:=\textsc{Zielonka}(V\setminus B)$ 
            \State $W_{1-\alpha } := W_{1-\alpha }\cup B$ 
        \EndIf
        \State \Return $W_{0},W_{1}$
    \EndFunction
\end{algorithmic}
```

# Solving with strategies

```{math}
\usepackage{algorithm}

\begin{algorithm}
\begin{algorithmic}[1]
\State {\textsc{SolveRecursive}}{$(V)$}
\State \textbf{if} {$ V = \emptyset$} \textbf{then} \Return {($\emptyset$,\ $\emptyset$,\ $\emptyset$,
 $\emptyset$)}
\State $m := \min(\{r(v) \mid v \in V\})$
\State $\alpha := m \mod 2$
\State $U := \{ v \in V \mid r(v) = m\}$
\State $A,S :=$ \textsc{Attr$_{\alpha}$}{$(U, V)$}
\State $W_{0}', W_{1}', S_{0}', S_{1}' := {\textsc{SolveRecursive}}{(V \setminus A)}$
\If {$W_{1-\alpha}' = \emptyset$} 
    \State $W_\alpha,W_{1 - \alpha} := A \cup W_{\alpha}', \emptyset$
    \State $S_\alpha := S \cup S_{\alpha}'\cup \{(u,v) \mid u \in U, v \in V\}$
    \State $S_{1 - \alpha} := \emptyset$
\Else
\State $B,S_B := $\textsc{Attr$_{1 - \alpha}$}{$(W_{1 - \alpha}', V)$}
\State $W_0, W_1, S_0, S_1 := {\textsc{SolveRecursive}}{(V \setminus B)}$
    \State $W_{1 - \alpha} := W_{1 - \alpha} \cup B$
    \State $S_{1 - \alpha} := S_{1 - \alpha}' \cup S_B \cup S_{1 - \alpha}$
    \State $\textcolor{red}{S_{\alpha} := S_{\alpha}' DELETE}$ \EndIf
\State \Return {$W_0, W_1, S_0, S_1$}
\end{algorithmic}
\caption{Solve recursive - with strategies}
\label{alg:ziel_symbolic_strategies}
\end{algorithm} 
```