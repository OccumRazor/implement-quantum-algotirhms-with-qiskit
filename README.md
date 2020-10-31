# implement-quantum-algotirhms-with-qiskit
Utilizing the python package, qiskit, published by IBM, (and of course other packages) to implement quantum algorithms.

This project is devoted to implement quantum algorithms, and currently the Grover search algorithm and the Deutsch-Jozsa algorithm have been finished. But note that the Grover search should be combined with a certain oracle while I just initialized the 'criteria quantum register' to a certain state, which can represent a binary number, and then compared it with the 'querry quantum register' quantum mechnically, so this can not strictly satisfy the description of the Grover search (searching an unstructured database), but what you need to do is simply exchange the corresponding code to a oracle_circuit.to_instruction() object and append it to the main circuit.

The later object of this project is to implement the Shor algorithm for integer factoring (other functions of this algorithm will not be given). Then, I'd like to try to implement quantum artificial intelligence algorithms, especially quantum neural networks.

If you viewed this project, I think it is reasonable to suppose that you know something about the fields of quantum computing, so the detailed introduction of each of the older algorithm will be omitted, but the papers have been cited inside some of the scripts. The furture works about quantum AI will contain detailed description since this is a active area of research and lots of the original papers mightbe quiet unfamiliar for some people.

Finally, I'll be quiet glad to share you my own papers and algorithms one day in the future!
