const redux = require('redux')
const createStore = redux.createStore;

const initState = {
  counter: 0
}

// Reducer : state will be updated here
const rootReducer = (state = initState, action) => {
  switch (action.type) {
    case 'INC_COUNTER':
      return { ...state, counter: state.counter + 1}
    default:
      return state
  }
}

// Store
const store = createStore(rootReducer);

// Subscription : will be informed on any state change
store.subscribe(() => {
  // will be called on any state change
  console.log('[Subscription', store.getState())
})

console.log(store.getState())

// Dispatching Action
store.dispatch({type: 'INC_COUNTER'})
console.log(store.getState())

