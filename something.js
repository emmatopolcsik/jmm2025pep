// name this Javascript file something.js
// execute by running this in a terminal:
//   node something.js
let a = 1, b = 1
console.log(a)
console.log(b)
Array.from({ length: 10 }, _ => {
  console.log(a+b)
  let _a = a
  a = b
  b = _a + b
})