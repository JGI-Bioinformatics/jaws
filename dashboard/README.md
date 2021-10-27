# This project (directory) was bootstrapped with [Create React App](https://github.com/facebook/create-react-app). 

### Setting up the development environment
As you change the react.js code, you can have a browser open (localhost:3000) and changes will automatically be updated there. 
You need to first follow these steps.

Do this on your **laptop** since you'll be using localhost.

1. Install node.js which will include npx and npm.  Here is an example for macOS.
https://nodejs.org/en/download/package-manager/#macos
note: don't use conda to install nodejs since this won't work

2. go into the dashboard project repository
git clone https://code.jgi.doe.gov/advanced-analysis/jaws.git
cd jaws/dashboard

3. source the secrets file (get this file from one of the JAWS react.js developers)
souce ~/.jaws-dashboard

4. install dependencies (these should not be stored in the repository since they are large files).
npm install

5. start the localhost server
npm start
```

# Some General Notes About npm Usage

## Available Scripts

In the project directory, you can run:

### `npm start`

Runs the app in the development mode.<br />
Open [https://jawsd.jgi.doe.gov:3003](https://jawsd.jgi.doe.gov:3003) to view it in the browser. The SSO server only works with **jgi.doe.gov** domain.

The page will reload if you make edits.<br />
You will also see any lint errors in the console.


The env.example defines certain default settings and env variables the app uses. You can customize these values by coping it to .env file, and edit the .env file accordingly. The the setting/env variables in .env file will be used by the app on a new dev server start.


### `npm test`
Launches the test runner in the interactive watch mode.<br />
See the section about [running tests](https://facebook.github.io/create-react-app/docs/running-tests) for more information.



### `npm run build`

Builds the app for production to the `build` folder.<br />
It correctly bundles React in production mode and optimizes the build for the best performance.

All the environmental variables the app needs (REACT_APP_*) need to be set here before running the build command.

The build is minified and the filenames include the hashes.<br />
Your app is ready to be deployed!

See the section about [deployment](https://facebook.github.io/create-react-app/docs/deployment) for more information.

### `npm run eject`

**Note: this is a one-way operation. Once you `eject`, you can’t go back!**

If you aren’t satisfied with the build tool and configuration choices, you can `eject` at any time. This command will remove the single build dependency from your project.

Instead, it will copy all the configuration files and the transitive dependencies (webpack, Babel, ESLint, etc) right into your project so you have full control over them. All of the commands except `eject` will still work, but they will point to the copied scripts so you can tweak them. At this point you’re on your own.

You don’t have to ever use `eject`. The curated feature set is suitable for small and middle deployments, and you shouldn’t feel obligated to use this feature. However we understand that this tool wouldn’t be useful if you couldn’t customize it when you are ready for it.

## Learn More

You can learn more in the [Create React App documentation](https://facebook.github.io/create-react-app/docs/getting-started).

To learn React, check out the [React documentation](https://reactjs.org/).

### Code Splitting

This section has moved here: https://facebook.github.io/create-react-app/docs/code-splitting

### Analyzing the Bundle Size

This section has moved here: https://facebook.github.io/create-react-app/docs/analyzing-the-bundle-size

### Making a Progressive Web App

This section has moved here: https://facebook.github.io/create-react-app/docs/making-a-progressive-web-app

### Advanced Configuration

This section has moved here: https://facebook.github.io/create-react-app/docs/advanced-configuration

### Deployment

This section has moved here: https://facebook.github.io/create-react-app/docs/deployment

### `npm run build` fails to minify

This section has moved here: https://facebook.github.io/create-react-app/docs/troubleshooting#npm-run-build-fails-to-minify
