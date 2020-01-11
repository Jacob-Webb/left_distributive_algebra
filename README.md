![lda image](references/img_on_3.png) <br />
# Left Distributive Algebras and Elementary Embeddings
I added a graph, node, and edge class from Qt to visually represent Dr. Scott Cramer's work on embeddings and left distributive algebras.  
### Motivation
During my time assisting Dr. Cramer I had the opportunity to study some set theory and learn more about large cardinal numbers and left-distributive algebras. Mostly I spent out time together understanding the code that he had written so that I could add the visual wrappers to it and to manipulate it in order to update it. For example, Dr. Cramer used arrays that sometimes went out of bounds so I adopted vectors. This was a simple change and didn't fundamentally alter any logic, but it allowed the program to keep from crashing at certain points. 
### Explanation
* Embedding: A mapping of one structure to another while maintaining first-order logic. So, given an embedding, j, and a structure N, j: N->M where M keeps without removing any first-order logic.
* Left-Distributive Algebra: If the left-(self)-distributivity is the property where a * (b * c) = (a * b) * (a * c), then given an elementary embedding j:Vλ → Vλ, Aj is the algebra of embeddings generated from j.
(Definition from Dr. Scott Cramer)
*  <a href="http://mathworld.wolfram.com/Rank.html" target="_blank">Rank</a>: The function from a set to an 
          <a href="http://mathworld.wolfram.com/OrdinalNumber.html" target="_blank">ordinal number</a>. The rank of a set is the smallest ordinal number that is greater 
          than all of the ranks of the elements of the set. 
          For our case an example of rank 3 would include: jjj, (jj)j, ... , j<sup>2</sup>j, j<sup>3</sup.  
          
### Tech Stack
* C++
* QT
### Resources
* <a href="references/Cramer_2019.pdf">Cramer S., Algebraic properties of elementary embeddings, Higher Recursion Theory and Set Theory, Singapore, 2019</a>
* <a href="references/Dougherty_Jech_97.pdf">Dougherty R., Jech T., Left-Distributive Embedding Algebras, Electronic Research Announcements of the American Mathematical Society, 1997, 3, 28-37
* <a href="references/Laver_Miller_2011.pdf">Laver R., Miller S.K., Left division in the free left distributive algebra on one generator, J. Pure Appl. Algebra, 2011,
215(3), 276–282</a>
* <a href="references/Laver_Miller_2013.pdf">Laver R., Miller S.K., The free one-generated left distributive algebra: Basics and a simplified proof of the division algorithm, Cent. Eur. J. Math, 2013, 11(12), 2150-2175</a>



