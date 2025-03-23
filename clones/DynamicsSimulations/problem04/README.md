---
title: "problem04"
author: "Vishal Paudel"
date: "2025/02/04"
---

> 4\. *Cross product: geometry vs components.*  
> 
>     a) The geometric definition of crossproduct is this $\vec{a} \times \vec{b}$ is a vector $\vec{c}$ with magnitude $\|\vec{a}\|\|\vec{b}\|\sin\theta_{ab}$ that is orthmorogonal to $\vec{a}$ and $vec{b}$ in the direction given by the right hand rule. Use this definition to findan alternative geometric definition involving projection (namely: project $\vec{b}$ onto theplane that is orthogonal to $\vec{a}$; then stretch it by $\vec{|a|}$; then rotate it $\pi/2$ around the $\vec{a}$ axis). Use that definition to show the distributive rule $\vec{a} \times (\vec{b} + \vec{c}) = \vec{a} \times \vec{b} + \vec{a} \times \vec{c}$.  
>     b) Then use the distributive rule to find the component formula for cross product,namely that $\vec{a} \times \vec{b} = (a_2 b_3 − a_3 b_2)\hat{e}_1 + (a_3 b_1 − a_1 b_3) \hat{e}_2 + (a_1 b_2 − a_2 b_1)\hat{e}_3$.  
> Note that, this distributive law implies that, for given $va \times \vec{v}$ is a linear operator.That is $va \times \vec{v}$ is a linear function fo $\vec{v}$. Later in the course we will use this to replacethe cross product with a tensor product. Hint: You can read about this in, say, theRuina/Pratap book (box 1.7).

The cross product $\vec{a} \times \vec{b}$ is defined geometrically as a vector $\vec{c}$ with magnitude $\|\vec{a}\|\|\vec{b}\|\sin \theta_{ab}$, direction orthogonal to both $\vec{a}$ and $\vec{b}$, and orientation following the right-hand rule. We can reformulate this using projection. First, project $\vec{b}$ onto the plane orthogonal to $\vec{a}$ using the formula $\vec{b}_{\perp} = \vec{b} - (\vec{b}\cdot\hat{a})\hat{a}$, where $\hat{a}$ is the unit vector in direction of $\vec{a}$. The magnitude of this projection is $\|\vec{b}\_{\perp}\| = \|\vec{b}\|\sin \theta$

Then stretch this vector by $\|\vec{a}\|$, giving $\|\vec{a}\|\vec{b}\_{\perp}$ with magnitude $\|\vec{a}\|\|\vec{b}\|\sin \theta$. Finally, rotate this vector by $\pi/2$ around the $\vec{a}$ axis, which preserves magnitude while making the result orthogonal to both $\vec{a}$ and $\vec{b}\_{\perp}$ in the right-hand rule direction. This construction yields a vector identical to the original definition.

To prove distributivity, $\vec{a} \times (\vec{b} + \vec{c}) = \vec{a} \times \vec{b} + \vec{a} \times \vec{c}$, we use this new geometric definition. When we project $(\vec{b} + \vec{c})$ onto the plane perpendicular to $\vec{a}$, linearity of projection gives us $(\vec{b} + \vec{c})\_{\perp} = \vec{b}\_{\perp} + \vec{c}\_{\perp}$. Stretching by $\|\vec{a}\|$ is also linear: $\|\vec{a}\|(\vec{b}\_{\perp} + \vec{c}\_{\perp}) = \|\vec{a}\|\vec{b}\_{\perp} + \|\vec{a}\|\vec{c}\_{\perp}$. Since rotation by $\pi/2$ is linear, $R\_{\pi/2}(\vec{b}\_{\perp} + \vec{c}\_{\perp}) = R\_{\pi/2}\vec{b}\_{\perp} + R\_{\pi/2}\vec{c}\_{\perp}$, proving the distributive property.

To derive the component formula, we use distributivity: $\vec{a} \times \vec{b} = (a\_1\hat{e}\_1 + a\_2\hat{e}\_2 + a\_3\hat{e}\_3) \times (b\_1\hat{e}\_1 + b\_2\hat{e}\_2 + b\_3\hat{e}\_3) = \sum\_{i,j} a\_ib\_j(\hat{e}\_i \times \hat{e}\_j)$. Using the standard basis cross products ($\hat{e}\_1 \times \hat{e}\_2 = \hat{e}\_3$, $\hat{e}\_2 \times \hat{e}\_3 = \hat{e}\_1$, $\hat{e}\_3 \times \hat{e}\_1 = \hat{e}\_2$, $\hat{e}\_i \times \hat{e}\_i = 0$, $\hat{e}\_j \times \hat{e}\_i = -(\hat{e}\_i \times \hat{e}\_j)$), we expand and group terms to get the final component formula: $\vec{a} \times \vec{b} = (a\_2b\_3 - a\_3b\_2)\hat{e}\_1 + (a\_3b\_1 - a\_1b\_3)\hat{e}\_2 + (a\_1b\_2 - a\_2b\_1)\hat{e}\_3$.
