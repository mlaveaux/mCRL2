A Telephone Book
----------------

**Contribution of this section**

#. use of functions with updates,
#. modelling considerations,
#. realistic verification process and its problems.

**New tools**:
none.

In this section we describe the model of a telephone book in
mCRL2.
We base our model on the following requirements:

#. A phone book shall store phone numbers.
#. It shall be possible to add and delete entries from a phone book.
#. It shall be possible to retrieve a phone number given a name.

By looking at these requirements, we identify the following entities:

* phone book
* phone number
* name

We start by giving abstract types for phone numbers and names; their concrete form
will be given later. For the phone book we decide that it is a mapping of names
to numbers. This gives us the following specification of sorts.

.. literalinclude:: files/phonebook1.mcrl2
   :language: mcrl2
   :lines: 4-6

As an mCRL2 user, you need to be aware that the types as given here, are
predetermined. This means that they cannot be extended on the fly. As a
consequence, all names and phone numbers that can ever be added to the phone
book must be known upfront.

If we again look at the requirements, our phone book must support the following
operations:

* ``addPhone`` adds a phone number for a name,
* ``delPhone`` deletes a phone number corresponding to a name, and
* ``findPhone`` finds the phone number corresponding to a name.

These operations will be the *actions* of our process. We need to decide on the
parameters that the actions are going to take. We assume that our process will
support a single phone book, i.e. the process itself models a phone book. It is
then natural to model ``addPhone`` with parameters ``Name`` and ``PhoneNumber``,
``delPhone`` with a ``Name``, and ``findPhone`` with a ``Name``. This gives rise
to the following action specification.

.. literalinclude:: files/phonebook1.mcrl2
   :language: mcrl2
   :lines: 12-15

We now need to take care that not every number is in every phone book. In order
to describe a phone book as a total function, we introduce a special phone number,
``p0``, to indicate that a name has no associated phone number.

.. literalinclude:: files/phonebook1.mcrl2
   :language: mcrl2
   :lines: 8-10

Given this decision, we can specify the empty phone book as the phone book that
maps every name to ``p0``.

.. code-block:: mcrl2

  lambda n: Name . p0;

In modelling the empty phone book we use lambda abstraction. In this expression
``lambda n: Name`` says that we are defining a function that takes arguments of
type ``Name``, and for each name produces ``p0`` as a result. As ``p0`` is of
type ``PhoneNumber``, ``lambda n: Name . p0`` describes a function of type
``Name ->PhoneNumber``, which is by definition equal to ``PhoneBook``.

Given a function ``b`` of type ``PhoneBook``, a name ``n`` and a phone number
``p``, we can set the value of ``n`` in ``b`` to ``p`` using the expression
``b[n -> p]``. This has as property that, for all names ``m != n``, ``b[n ->
p](m) = b(m)``, and ``b[n -> p](n) = p``.

Using the above ingredients, we can model a simple phone book using the
following specification, that is also available as
:download:`phonebook1.mcrl2 <files/phonebook1.mcrl2>`.

.. literalinclude:: files/phonebook1.mcrl2
   :language: mcrl2

.. admonition:: Exercise

  There are some obvious flaws in the phone book that we have specified. Can you
  find and explain them?

.. admonition:: Hint
   :class: collapse

   Think about the special meaning of phone number ``p0``, and explain the
   ``findPhone`` function to yourself.

.. admonition:: Solution
   :class: collapse

   In this specification, the "special" phone number ``p0`` can be assigned to a
   name freely. Furthermore, a ``findPhone`` action can be performed, but the
   actual phone number is never reported.

.. admonition:: Exercise

   Fix these issues in the above specification.

.. admonition:: Solution
   :class: collapse

   Preventing the assignment of ``p0`` to a name can easily prevented by
   guarding the ``addPhone`` action with ``p != p0``. Fixing the second issue
   requires some more thought. There are two possible ways around fixing the
   problem. We can either assume that reporting of the result is immediate, and
   add the resulting phone number as a parameter to the ``findPhone`` action, or
   we can assume that querying for a phone number is asynchronous, and my take
   time, and split the query into the action initiating the query
   (``findPhone``) and an action reporting the result, e.g. ``reportPhone``.

   The first approach is suitable when, e.g., modelling a phone book that is a
   library in a synchronous program in, say, C or Java. In that case indeed the
   program pointer of the calling program will not change before the result has
   been returned, making a model in which reporting the result a faithful
   representation of reality. This variation is implemented in the following
   specification (available as :download:`phonebook2a.mcrl2
   <files/phonebook2a.mcrl2>`). The result is obtained instantaneously. Watch
   the extra parameter to ``addPhone``.

   .. literalinclude:: files/phonebook2a.mcrl2
      :language: mcrl2

   If we are, e.g., modelling a phone book that is a web service, where a client
   performs a request, and in the meantime may do other kinds of actions like
   sending requests to other web services, the previous approach provides too
   coarse an abstraction. In this case it is more accurate to use the second
   approach, in which performing the query and obtaining the result are truely
   decoupled. This variation is given by the following specification (available
   as :download:`phonebook2b.mcrl2 <files/phonebook2b.mcrl2>`).

   .. literalinclude:: files/phonebook2b.mcrl2
      :language: mcrl2

In the rest of this tutorial we will stick to the specification with
asynchronous reporting (the second variation in the previous exercise).

In complex specifications, it can be convenient to introduce additional
functions, with descriptive names, that take care of the modifications of
parameters that is done in the process. As a bonus this usually makes it easier
to change the data structures used in a specification.

.. admonition:: Exercise

   Modify the specification in :download:`phonebook2b.mcrl2
   <files/phonebook2b.mcrl2>` by adding functions ``emptybook, add_phone,
   del_phone`` and ``find_phone`` with the following signatures.

   .. code-block:: mcrl2

      map  emptybook: PhoneBook;
           add_phone: PhoneBook # Name # PhoneNumber -> PhoneBook;
           del_phone: PhoneBook # Name -> PhoneBook;
           find_phone: PhoneBook # Name -> PhoneNumber;

   Define equations implementing the above operations.

.. admonition:: Solution
   :class: collapse

   A solution to the above exercise is given by the following specification
   (also available as :download:`phonebook3.mcrl2 <files/phonebook3.mcrl2>`).

   .. literalinclude:: files/phonebook3.mcrl2
      :language: mcrl2

It should be noted that, instead of using a function of names to phone numbers,
we could also have modelled the phone book using a set of pairs of names and
phone numbers. A model using sets is likely to become complicated in this case.

.. admonition:: Exercise

   Modify the specification in :download:`phonebook3.mcrl2
   <files/phonebook3.mcrl2>` such that it uses a set of pairs of names and phone
   numbers instead of function from names to phone numbers to store the phone
   numbers internally.

.. admonition:: Solution
   :class: collapse

   A possible solution to this exercise is given in
   :download:`phonebook4.mcrl2 <files/phonebook4.mcrl2>`. Note that the function
   ``find_phone`` cannot be implemented using sets, because no concrete elements
   can be taken from the set. Therefore, the functionality of ``find_phone`` is
   modelled using the ``sum`` operator on a process level.

In the rest of this chapter we stick to the model in which functions occur
directly in the specification (:download:`phonebook2b.mcrl2
<files/phonebook2b.mcrl2>`). We are going to check whether our model makes
sense using the µ-calculus.

A suitable property for our specification is the following:

  "if a name ``n`` with phone number ``p`` is added to the phone book, and a
  lookup of name ``n`` is performed, then phone number ``p`` should be reported,
  provided that in the meantime no other phone number has been added for name
  ``n``, and the phone number for name ``n`` has not been deleted".

.. admonition:: Exercise

   Write the above propery in the µ-calculus.

.. admonition:: Solution
   :class: collapse

   A solution (available as :download:`phonebook1.mcf <files/phonebook1.mcf>`)
   is the following:

   .. literalinclude:: files/phonebook1.mcf
      :language: mcrl2

.. admonition:: Exercise

   Try to verify this property using the tools. What do you observe?

.. admonition:: Hint
   :class: collapse

   You need to use the tools :ref:`tool-mcrl22lps`, :ref:`tool-lps2pbes` and
   :ref:`tool-pbessolve`.

.. admonition:: Solution
   :class: collapse

   We try checking this property using the following command::

     $ mcrl22lps phonebook2b.mcrl2 | lps2pbes -f phonebook1.mcf | pbessolve
     [20:54:11.222 error]   Cannot find a term of sort Name

   Observe that the tools fail to verify this requirement because of an error.

   The tool is telling you that it wants to find some representative term of
   sort ``Name``, but is not able to do so. This indeed makes sense as we have
   not given a specification of names and phone numbers yet.

.. admonition:: Exercise

   Fix the specification by making sure there are 3 names, and 3 phone numbers
   (apart from the special phone number ``p0``). Rerun the verification on the
   fixed specification.

.. admonition:: Solution
   :class: collapse

   The following is a fixed specification (also available as
   :download:`phonebook5.mcrl2 <files/phonebook5.mcrl2>`).

   .. literalinclude:: files/phonebook5.mcrl2
      :language: mcrl2

   The specification is now easily checked using the following sequence of
   commands::

     $ mcrl22lps phonebook5.mcrl2 | lps2pbes -f phonebook1.mcf | pbessolve
     true

.. admonition:: Exercise

   Verify whether the following property holds for :download:`phonebook5.mcrl2
   <files/phonebook5.mcrl2>`.

     "if a name ``n`` with phone number ``p`` is added to the phone book, and a
     lookup of name ``n`` is performed, then phone number ``p`` should be
     reported, provided that in the meantime the phone number for name ``n`` has
     not been deleted".

   You first need to formalise this property as a  µ-calculus formula, and then
   verify whether it holds. Explain the outcome of the verification.

.. admonition:: Solution
   :class: collapse

   The following formula (available as :download:`phonebook2.mcf <files/phonebook2.mcf>`)
   formalises this property.

   .. literalinclude:: files/phonebook2.mcf
      :language: mcrl2

   It can be verified using the commands::

     $ mcrl22lps phonebook5.mcrl2 | lps2pbes -f phonebook2.mcf | pbessolve
     false

   Observe that this verification fails because ``addPhone`` allows you to add
   a phone number for a person that already has a phone number. If a new phone
   number is added for such a person, the original phone number is overwritten.

.. admonition:: Exercise

   Modify the specification in :download:`phonebook5.mcrl2
   <files/phonebook5.mcrl2>` such that ``addPhone(n,p)`` can only be executed if
   no phone number for name ``n`` is known. Furthermore, extend the
   specification with and action ``changePhone`` with the following signature

   .. code-block:: mcrl2

      changePhone: Name # PhoneNumber

   such that ``changePhone(n,p)`` can only be executed if ``n`` already has a
   phone number, and that afterwards the phone number of ``n`` has been updated
   to ``p``. Save the result as :file:`phonebook6.mcrl2`

.. admonition:: Solution
   :class: collapse

   A sample solution is the following (:download:`phonebook6.mcrl2
   <files/phonebook6.mcrl2>`)

   .. literalinclude:: files/phonebook6.mcrl2
      :language: mcrl2

.. admonition:: Exercise

   Verify whether your new specification satisfies the property you formulated
   before. Explain the outcome.

.. admonition:: Solution

   For our version of the property, the verification delivers the following
   result::

     $ mcrl22lps phonebook6.mcrl2 | lps2pbes -f phonebook2.mcf | pbessolve
     false

.. admonition:: Exercise

   If the verification in the previous exercise failed, think about the influence
   of the ``changePhone`` action on the validity of the property you are trying
   to check. Describe the changed property in natural language, give the modal
   µ-calculus formula, and do the verification.

.. admonition:: Solution

   Our modified property is the following (:download:`phonebook3.mcf <files/phonebook3.mcf>`).

   .. literalinclude:: files/phonebook3.mcf
      :language: mcrl2

   Now the following verification succeeds::

     $ mcrl22lps phonebook6.mcrl2 | lps2pbes -f phonebook3.mcf | pbessolve
     true

As an alternative approach to the verification, you can make the phonebook that
is used internally available externally through a ``getPhoneBook`` action.
A specification in which this modification has been made is the following
(:download:`phonebook7.mcrl2 <files/phonebook7.mcrl2>`).

.. literalinclude:: files/phonebook7.mcrl2
   :language: mcrl2

Our last property (:download:`phonebook3.mcf <files/phonebook3.mcf>`) can now
also be formulated as follows (:download:`phonebook4.mcf <files/phonebook4.mcf>`)
using the internally stored phonebook.

.. literalinclude:: files/phonebook4.mcf
   :language: mcrl2

The verification result is still the same::

  $ mcrl22lps phonebook7.mcrl2 | lps2pbes -f phonebook4.mcf | pbessolve
  true
